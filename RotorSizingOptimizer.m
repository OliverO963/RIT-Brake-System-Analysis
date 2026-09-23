%% ================================================================
%  Rotor Sizing Optimizer - target mass & convective area per rotor
% ================================================================
%  Consumes:  a drive-cycle CSV (either force-based, e.g.
%             Regenerative_Braking_Forces_8_Laps.csv, or the legacy
%             torque-based MotorTorque.csv format - auto-detected)
%             VehicleParameters.xlsx  (the single shared parameter
%             sheet also read by loadVehicleParams.m; all rotor-sizing
%             constants live on its Parameters tab, read BY NAME)
%
%  Produces:  m_target and A_target for the FRONT and REAR rotor,
%             an existence-proof geometry, a mass/area trade-off
%             table, the achievable geometric envelope, the active
%             limits, and an h_w / T_target sensitivity study.
%             Written to RotorSizingResults.xlsx for the
%             generative-design tool.
%
%  The thermal physics is run_sim_opt() from BrakeCoeffOptimizer.m,
%  reused unchanged (see run_sim_rotor below). The lsqnonlin fit is
%  NOT re-run: the fitted h_w and PadFrac coefficients are INPUTS.
%
%  FORCE-BASED FORMAT CORRECTION: FrictionBrakeForce_N,
%  RequiredBrakingMagnitude_N, AvailableRegenMagnitude_N and
%  AppliedRegenForce_N are longitudinal forces at the contact patch
%  (not pedal input, not torque). A braking event is active whenever
%  RequiredBrakingMagnitude_N > 0. FrictionBrakeForce_N is used
%  directly as the friction-heat source (see ingest_new_format).
%
%  Requires: base MATLAB only (no toolboxes).
% ================================================================

clc; clear; close all

%% ================== USER INPUTS ==================
cfg.paramFile   = 'VehicleParameters.xlsx';
cfg.cycleFile   = 'Regenerative_Braking_Forces_8_Laps.csv';
cfg.resultsFile = 'RotorSizingResults.xlsx';
cfg.makePlots   = true;
cfg.verbose     = true;

P = load_rotor_params(cfg.paramFile);
sides = {'front','rear'};

%% ================== LOAD AND VALIDATE THE DRIVE CYCLE ==================
C = readtable(cfg.cycleFile);
cols = C.Properties.VariableNames;
isNewFormat = all(ismember({'Lap','Time_s','Speed_mps','RequiredBrakingMagnitude_N', ...
    'AvailableRegenMagnitude_N','AppliedRegenForce_N','FrictionBrakeForce_N', ...
    'AeroDragMagnitude_N','LongitudinalForce_N'}, cols));
isOldFormat = all(ismember({'time','dt_s','speed_mps','motor_rpm', ...
    'drive_torque_per_motor_Nm','regen_torque_per_motor_Nm', ...
    'torque_per_motor_Nm','total_motor_torque_Nm','friction_brake_force_N'}, cols));
if ~isNewFormat && ~isOldFormat
    error('RotorSizing:BadCycle', ...
        '%s matches neither the force-based format nor the legacy torque-based format.', ...
        cfg.cycleFile);
end

if isNewFormat
    [CYC, D, V, n] = ingest_new_format(C, P);
else
    [CYC, D, V, n] = ingest_old_format(C, P);
end
lapEdges = [find(CYC.lapStart); n+1];

if cfg.verbose
    print_diagnostics(V, D, P, cfg, sides, isNewFormat);
end

%% ================== ROTOR GEOMETRY DEFINITIONS ==================
G.front = make_geom(P.Do_front_mm/1000, P.Di_front_mm/1000, P.t_min_front_mm/1000, ...
    P.t_max_front_mm/1000, P.s_min_front_mm/1000, P.phi_max_front, P.cut_type_front, ...
    P.include_edge_area, P.rho_rotor, P.edge_margin_mm/1000, P.enforce_pattern_packing);
G.rear  = make_geom(P.Do_rear_mm/1000,  P.Di_rear_mm/1000,  P.t_min_rear_mm/1000, ...
    P.t_max_rear_mm/1000,  P.s_min_rear_mm/1000,  P.phi_max_rear,  P.cut_type_rear, ...
    P.include_edge_area, P.rho_rotor, P.edge_margin_mm/1000, P.enforce_pattern_packing);

hw.front = [P.hw_x1_front, P.hw_b1_front];
hw.rear  = [P.hw_x1_rear,  P.hw_b1_rear];
TambK = P.TambC_fallback + 273.15;
TinitK = (P.T_init_F - 32)*(5/9) + 273.15;
padp  = [P.padfrac_x2, P.padfrac_b2];

%% ================== SOLVE EACH ROTOR ==================
RES = struct();
for s = 1:2
    sd = sides{s};
    g  = G.(sd);
    sim = @(m,A,hs) run_sim_rotor(CYC, sd, m, A, hw.(sd)*1, padp, P, TambK, TinitK, hs);
    pk  = @(m,A,hs) max(K2F(sim(m,A,hs)));

    mlo = g.rho*g.tmin*(1-g.phimax)*g.Aring;
    mhi = g.rho*g.tmax*g.Aring;
    mgrid = linspace(mlo, mhi, P.n_mass_grid)';

    Astar = zeros(size(mgrid)); Agmax = zeros(size(mgrid)); Agmin = zeros(size(mgrid));
    tAtMax = zeros(size(mgrid)); ArAtMax = zeros(size(mgrid)); LAtMax = zeros(size(mgrid));
    pkMax  = zeros(size(mgrid)); pkMin = zeros(size(mgrid)); reach = false(size(mgrid));
    for k = 1:numel(mgrid)
        Astar(k) = area_frontier(pk, mgrid(k), P.T_target_F, P.A_bracket_lo, ...
                                 P.A_bracket_hi, P.area_bisect_iters, 1);
        [Agmax(k), tAtMax(k), ArAtMax(k), LAtMax(k)] = geo_area_max(g, mgrid(k));
        Agmin(k) = geo_area_min(g, mgrid(k));
        pkMax(k) = pk(mgrid(k), Agmax(k), 1);
        pkMin(k) = pk(mgrid(k), Agmin(k), 1);
        reach(k) = Agmax(k) >= Astar(k);
    end

    kOpt = find(reach, 1, 'first');
    if isempty(kOpt)
        [~, kBest] = min(pkMax);
        warning('RotorSizing:Infeasible', ...
            ['%s: T_target = %.1f degF is NOT achievable anywhere in the geometric envelope. ' ...
             'Lowest achievable peak = %.1f degF (at m = %.4f kg, A = %.5f m^2).'], ...
            upper(sd), P.T_target_F, min(pkMax), mgrid(kBest), Agmax(kBest));
        kOpt = kBest;
        RES.(sd).feasible = false;
    else
        RES.(sd).feasible = true;
    end

    m_t = mgrid(kOpt);
    A_t = Astar(kOpt);
    [Ag, tg, Arg, Lg] = geo_area_max(g, m_t);

    % --- Existence proof: a concrete (t, Ar, L) that is BOTH >= A_t and
    %     buildable as a single row of drilled holes / radial slots.
    proof = existence_proof(g, m_t, A_t);

    % --- Re-simulate the recommendation ---
    Tk = sim(m_t, max(A_t, proof.A), 1);
    TF = K2F(Tk);
    lapPeak = zeros(numel(lapEdges)-1,1);
    for k = 1:numel(lapEdges)-1
        lapPeak(k) = max(TF(lapEdges(k):lapEdges(k+1)-1));
    end
    TkAt = sim(m_t, A_t, 1);

    % --- Explicit-Euler stability guard ---
    cpMin = (0.0005*TinitK + 0.2813)*1000;
    hwMax = hw.(sd)(1)*max(velx) + hw.(sd)(2);
    euler = max(CYC.tstep) * hwMax * max(Ag, proof.A) / (m_t * cpMin);
    if euler > P.euler_warn_threshold
        warning('RotorSizing:Euler', ...
            '%s: explicit-Euler number dt*h*A/(m*cp) = %.3f (warn > %.2f). Sub-step the cycle.', ...
            upper(sd), euler, P.euler_warn_threshold);
    end

    % --- Active limits ---
    lim = {};
    if abs(tg - g.tmin) < 1e-9,                       lim{end+1} = 't_min';   end %#ok<AGROW>
    if abs(tg - g.tmax) < 1e-9,                       lim{end+1} = 't_max';   end %#ok<AGROW>
    if abs(Arg - g.phimax*g.Aring)/g.Aring < 1e-6
        if g.packingLimited
            lim{end+1} = sprintf('phi_max (PACKING-LIMITED to %.3f, spec %.3f)', ...
                                 g.phimax_pack, g.phimax_spec); %#ok<AGROW>
        else
            lim{end+1} = 'phi_max'; %#ok<AGROW>
        end
    end
    if Lg >= g.cutC*Arg/g.smin - 1e-9 && Arg > 0,     lim{end+1} = 's_min';   end %#ok<AGROW>
    if RES.(sd).feasible && kOpt > 1,                 lim{end+1} = 'T_target'; end %#ok<AGROW>
    if RES.(sd).feasible && kOpt == 1
        lim{end+1} = 'T_target SLACK (mass floor governs)'; %#ok<AGROW>
    end

    % --- Hole-vs-2t area-gain criterion ---
    % dA/dAr at fixed t with L = cutC*Ar/s_min is (cutC*t/s_min - 2).
    gainPossible = (g.cutC*g.tmax/g.smin - 2) > 0;
    tBreakEven   = 2*g.smin/g.cutC;

    RES.(sd).g = g;  RES.(sd).mgrid = mgrid; RES.(sd).Astar = Astar;
    RES.(sd).Agmax = Agmax; RES.(sd).Agmin = Agmin; RES.(sd).reach = reach;
    RES.(sd).tAtMax = tAtMax; RES.(sd).ArAtMax = ArAtMax; RES.(sd).LAtMax = LAtMax;
    RES.(sd).pkMax = pkMax; RES.(sd).pkMin = pkMin;
    RES.(sd).m_target = m_t; RES.(sd).A_target = A_t;
    RES.(sd).A_geo_max_at_target = Ag; RES.(sd).A_geo_min_at_target = geo_area_min(g,m_t);
    RES.(sd).proof = proof; RES.(sd).TF = TF; RES.(sd).lapPeak = lapPeak;
    RES.(sd).peakF = max(TF); RES.(sd).endF = TF(end);
    RES.(sd).peakF_atAtarget = max(K2F(TkAt));
    RES.(sd).m_range = [mlo mhi]; RES.(sd).activeLimits = strjoin(lim, ', ');
    RES.(sd).euler = euler;
    RES.(sd).gainPossible = gainPossible; RES.(sd).tBreakEven = tBreakEven;

    % --- Sensitivity ---
    sens = struct('case',{},'m_target',{},'A_target',{});
    cases = {'baseline',1,P.T_target_F; 'h_w -20%',0.8,P.T_target_F; 'h_w +20%',1.2,P.T_target_F; ...
             'T_target -25F',1,P.T_target_F-25; 'T_target +25F',1,P.T_target_F+25};
    for q = 1:size(cases,1)
        hs = cases{q,2}; Tt = cases{q,3};
        kk = [];
        for k = 1:numel(mgrid)
            a = area_frontier(pk, mgrid(k), Tt, P.A_bracket_lo, P.A_bracket_hi, P.area_bisect_iters, hs);
            if Agmax(k) >= a, kk = k; break; end
        end
        if isempty(kk), sens(q).case = cases{q,1}; sens(q).m_target = NaN; sens(q).A_target = NaN;
        else
            sens(q).case = cases{q,1}; sens(q).m_target = mgrid(kk);
            sens(q).A_target = area_frontier(pk, mgrid(kk), Tt, P.A_bracket_lo, ...
                                             P.A_bracket_hi, P.area_bisect_iters, hs);
        end
    end
    RES.(sd).sens = sens;
end

%% ================== REPORT ==================
fprintf('\n================ ROTOR TARGETS (T_target = %.1f degF) ================\n', P.T_target_F);
for s = 1:2
    sd = sides{s}; R = RES.(sd); g = R.g; pr = R.proof;
    fprintf('\n--- %s rotor ---\n', upper(sd));
    fprintf('A_ring (fixed annulus face) : %.6f m^2   radial width %.2f mm\n', ...
        g.Aring, (g.Do-g.Di)/2*1000);
    fprintf('achievable mass envelope    : %.4f - %.4f kg\n', R.m_range(1), R.m_range(2));
    if g.packingLimited
        fprintf('*** phi_max reduced %.3f -> %.3f: no single-row %s pattern with a %.2f mm radial\n', ...
            g.phimax_spec, g.phimax_pack, g.cutKind, g.edge*1000);
        fprintf('    margin removes more than %.1f%% of the face at this radial width (%.2f mm).\n', ...
            100*g.phimax_pack, (g.Do-g.Di)/2*1000);
    end
    fprintf('m_target                    : %.4f kg\n', R.m_target);
    fprintf('A_target = A*(m_target)     : %.5f m^2   (thermal minimum)\n', R.A_target);
    fprintf('geometric area at m_target  : %.5f - %.5f m^2\n', R.A_geo_min_at_target, R.A_geo_max_at_target);
    fprintf('EXISTENCE PROOF (feasibility only - NOT the design):\n');
    fprintf('   t = %.3f mm, Ar = %.6f m^2 (phi = %.3f), L = %.4f m\n', ...
        pr.t*1000, pr.Ar, pr.Ar/g.Aring, pr.L);
    fprintf('   pattern: %d x %.2f mm %s, single row on a %.2f mm pitch circle\n', ...
        pr.N, pr.feat*1000, pr.kind, pr.Dm*1000);
    fprintf('   recomputed m = %.4f kg, A = %.5f m^2\n', pr.m_check, pr.A);
    fprintf('resulting temps at the proof geometry: peak %.1f degF, end-of-run %.1f degF\n', ...
        R.peakF, R.endF);
    fprintf('   per-lap peaks: %.1f - %.1f degF\n', min(R.lapPeak), max(R.lapPeak));
    fprintf('peak at exactly A_target    : %.1f degF (target %.1f)\n', R.peakF_atAtarget, P.T_target_F);
    fprintf('active limits               : %s\n', R.activeLimits);
    fprintf('explicit-Euler number       : %.4f\n', R.euler);
    fprintf('cutting raises area?        : %s (break-even thickness %.2f mm for this cut type)\n', ...
        ternary(R.gainPossible,'yes','no'), R.tBreakEven*1000);
    fprintf('sensitivity:\n');
    for q = 1:numel(R.sens)
        fprintf('   %-14s m_target = %.4f kg, A_target = %.5f m^2\n', ...
            R.sens(q).case, R.sens(q).m_target, R.sens(q).A_target);
    end
end

%% ================== WRITE RESULTS FILE ==================
Summary = table();
for s = 1:2
    sd = sides{s}; R = RES.(sd);
    row = table(string(sd), R.m_target, R.A_target, R.A_geo_min_at_target, R.A_geo_max_at_target, ...
        R.proof.t*1000, R.proof.Ar, R.proof.L, R.proof.N, R.proof.feat*1000, R.proof.A, ...
        R.peakF, R.endF, P.T_target_F, R.m_range(1), R.m_range(2), string(R.activeLimits), R.feasible, ...
        'VariableNames', {'Rotor','m_target_kg','A_target_m2','A_geo_min_m2','A_geo_max_m2', ...
        'proof_t_mm','proof_Ar_m2','proof_L_m','proof_N_cuts','proof_feature_mm','proof_A_m2', ...
        'peak_degF','end_degF','T_target_degF','m_env_min_kg','m_env_max_kg','active_limits','temp_feasible'});
    Summary = [Summary; row]; %#ok<AGROW>
end
writetable(Summary, cfg.resultsFile, 'Sheet', 'Targets');

for s = 1:2
    sd = sides{s}; R = RES.(sd);
    Tr = table(R.mgrid, R.Astar, R.Agmin, R.Agmax, R.tAtMax*1000, R.ArAtMax/R.g.Aring, ...
        R.LAtMax, R.reach, R.pkMax, R.pkMin, 'VariableNames', ...
        {'m_kg','A_required_m2','A_geo_min_m2','A_geo_max_m2','t_at_Amax_mm', ...
         'phi_at_Amax','L_at_Amax_m','geometrically_reachable','peak_at_A_geo_max_degF', ...
         'peak_at_A_geo_min_degF'});
    writetable(Tr, cfg.resultsFile, 'Sheet', ['Tradeoff_' sd]);
    Sn = struct2table(R.sens);
    writetable(Sn, cfg.resultsFile, 'Sheet', ['Sensitivity_' sd]);
end

Diag = table(fieldnames(rmfield(D,{'perLap_Efric_front_kJ','perLap_Efric_rear_kJ'})), ...
    struct2cell(rmfield(D,{'perLap_Efric_front_kJ','perLap_Efric_rear_kJ'})), ...
    'VariableNames', {'Metric','Value'});
writetable(Diag, cfg.resultsFile, 'Sheet', 'Diagnostics');
writetable(table((1:numel(D.perLap_Efric_front_kJ))', D.perLap_Efric_front_kJ, D.perLap_Efric_rear_kJ, ...
    'VariableNames', {'Lap','Efric_front_kJ','Efric_rear_kJ'}), cfg.resultsFile, 'Sheet', 'PerLap');
fprintf('\nResults written to %s\n', cfg.resultsFile);

%% ================== PLOTS ==================
if cfg.makePlots
    figure('Name','Rotor temperature at recommended targets');
    for s = 1:2
        sd = sides{s}; R = RES.(sd);
        subplot(2,1,s);
        plot(CYC.t, R.TF, 'b-'); hold on;
        yline(P.T_target_F, 'r--', 'T\_target');
        xlabel('Time (s)'); ylabel('Rotor temp (degF)'); grid on;
        title(sprintf('%s rotor: m = %.4f kg, A = %.5f m^2 (peak %.1f degF)', ...
            sd, R.m_target, max(R.A_target,R.proof.A), R.peakF), 'Interpreter','none');
    end

    for s = 1:2
        sd = sides{s}; R = RES.(sd);
        figure('Name',sprintf('Mass-area plane: %s', sd));
        fill([R.mgrid; flipud(R.mgrid)], [R.Astar; flipud(R.Agmax)], [0.85 0.95 0.85], ...
            'EdgeColor','none','DisplayName','feasible region'); hold on;
        plot(R.mgrid, R.Astar, 'r-',  'LineWidth', 2, 'DisplayName','temperature frontier A*(m)');
        plot(R.mgrid, R.Agmax, 'b-',  'LineWidth', 2, 'DisplayName','geometric envelope A_{geo,max}(m)');
        plot(R.mgrid, R.Agmin, 'b--', 'LineWidth', 1, 'DisplayName','geometric floor A_{geo,min}(m)');
        plot(R.m_target, max(R.A_target,R.proof.A), 'kp', 'MarkerSize', 14, ...
            'MarkerFaceColor','y', 'DisplayName','recommended point');
        xlabel('Rotor mass (kg)'); ylabel('Convective area (m^2)'); grid on;
        legend('Location','best'); title(sprintf('%s rotor mass-area plane', sd), 'Interpreter','none');
    end
end

%% ================================================================
%  LOCAL FUNCTIONS
% ================================================================


function [CYC, D, V, n] = ingest_new_format(C, P)
%INGEST_NEW_FORMAT Force-based drive-cycle format (Lap, Time_s, AccelX_g,
% MotorRPM, Speed_mps, LongitudinalForce_N, AeroDragMagnitude_N,
% RequiredBrakingMagnitude_N, AvailableRegenMagnitude_N,
% AppliedRegenForce_N, FrictionBrakeForce_N).
%
% CORRECTION APPLIED: FrictionBrakeForce_N, RequiredBrakingMagnitude_N,
% AvailableRegenMagnitude_N and AppliedRegenForce_N are all longitudinal
% forces at the CONTACT PATCH (whole vehicle), not pedal input and not
% torques. A braking event is active whenever RequiredBrakingMagnitude_N
% is nonzero (it is clipped at 0, so "nonzero" means "> 0").
%
% Because FrictionBrakeForce_N is the ACTUAL friction force already net of
% regen (verified: FrictionBrakeForce_N + AppliedRegenForce_N =
% -RequiredBrakingMagnitude_N to within 1e-6 N on the supplied file), the
% friction heat at each corner is now a DIRECT measurement -
% |FrictionBrakeForce_N| * distance, split by Tbias_front and the 0.5 L/R
% split - rather than an estimate backed out of a KE-loss/regen
% subtraction. This removes the old max(Ecorner-regen,0) clamp entirely
% (there is nothing left to clamp) and no longer needs the regen
% front/rear split for the heat calculation itself.
n = height(C);
sides = {'front','rear'};

V = struct();
V.rows       = n;
V.nan_count  = sum(sum(ismissing(C)));
t     = C.Time_s;
tstep = [0; diff(t)];
V.span         = t(end) - t(1);
V.dt_median    = median(tstep(2:end));
V.dt_uniform   = all(abs(tstep(2:end) - V.dt_median) < 1e-9);

velx = C.Speed_mps;
velx(abs(velx) > P.velx_threshold) = 0;
velx(isnan(velx)) = 0;
velx(velx < 0)    = 0;
V.n_speed_cleaned = sum(velx ~= C.Speed_mps);

lapStart = [true; diff(C.Lap) ~= 0];
V.lap_count = sum(lapStart);
lapRows     = diff([find(lapStart); n+1]);
V.lap_rows_equal = all(lapRows == lapRows(1));
V.lap_rows_first = lapRows(1);

% ---- Wheel-radius consistency check (diagnostic only) ------------------
omega_wheel = C.MotorRPM / P.gear_ratio * 2*pi/60;
r_implied   = velx ./ max(omega_wheel, eps);
V.r_implied_mean = mean(r_implied);
V.r_implied_vs_WheelR_pct = 100*(V.r_implied_mean/P.WheelR - 1);

RB = C.RequiredBrakingMagnitude_N;
FB = C.FrictionBrakeForce_N;
AR = C.AppliedRegenForce_N;
AV = C.AvailableRegenMagnitude_N;
LF = C.LongitudinalForce_N;
AD = C.AeroDragMagnitude_N;

active = RB > 0;
V.active_steps = sum(active);

% ---- Identity checks on the corrected columns (cheap, high-value) ------
V.identity_FB_AR_RB_maxerr = max(abs(FB + AR + RB));
V.identity_LF_RB_AD_maxerr_braking = max(abs(LF(active) - (-RB(active) - AD(active))));
V.n_active_with_AccelX_nonneg = sum(active & (C.AccelX_g >= 0));

dstep = zeros(n,1);
dstep(2:end) = (velx(1:end-1) + velx(2:end))/2 .* tstep(2:end);
Efric_total = abs(FB) .* dstep;      % vehicle-total friction energy, this step

Tbias.front = P.Tbias_front;
Tbias.rear  = 1 - P.Tbias_front;
for s = 1:2
    sd = sides{s};
    CYC.(sd).Efric = Efric_total * Tbias.(sd) * 0.5;   % front/rear x L/R split
end
CYC.tstep    = tstep;
CYC.velx     = velx;
CYC.active   = active;
CYC.t        = t;
CYC.lapStart = lapStart;

% ---- Diagnostics ----
D.rows = n; D.active_steps = V.active_steps; D.laps = V.lap_count;
D.regen_applied_total_kJ   = sum(abs(AR).*dstep)/1e3;
D.regen_available_total_kJ = sum(AV.*dstep)/1e3;
D.required_braking_total_kJ = sum(RB.*dstep)/1e3;
for s = 1:2
    sd = sides{s};
    D.(['Efric_' sd '_kJ']) = sum(CYC.(sd).Efric(active))/1e3;
end
D.friction_share_of_required_pct = 100*sum(Efric_total(active)) / max(D.required_braking_total_kJ*1e3,1);
D.identity_FB_AR_RB_maxerr = V.identity_FB_AR_RB_maxerr;

lapEdgesLocal = [find(lapStart); n+1];
D.perLap_Efric_front_kJ = zeros(numel(lapEdgesLocal)-1,1);
D.perLap_Efric_rear_kJ  = zeros(numel(lapEdgesLocal)-1,1);
for k = 1:numel(lapEdgesLocal)-1
    rows = lapEdgesLocal(k):lapEdgesLocal(k+1)-1;
    D.perLap_Efric_front_kJ(k) = sum(CYC.front.Efric(rows) .* active(rows))/1e3;
    D.perLap_Efric_rear_kJ(k)  = sum(CYC.rear.Efric(rows)  .* active(rows))/1e3;
end
end


function [CYC, D, V, n] = ingest_old_format(C, P)
%INGEST_OLD_FORMAT Legacy torque-based drive-cycle format (MotorTorque.csv-
% style: time, dt_s, speed_mps, motor_rpm, drive/regen/total torque per
% motor, friction_brake_force_N as PEDAL FORCE only). Friction heat is
% estimated by subtracting a per-corner regen split from the KE-loss chain,
% clamped at zero. Kept for backward compatibility with older drive cycles;
% new work should use the force-based format (see ingest_new_format).
n = height(C);

V = struct();
V.rows        = n;
V.nan_count   = sum(sum(ismissing(C)));
V.t_raw       = C.time;
V.dt_s        = C.dt_s;
V.span_raw    = C.time(end) - C.time(1);
V.x4_maxerr   = max(abs(C.total_motor_torque_Nm - 4*C.torque_per_motor_Nm));
V.sum_maxerr  = max(abs(C.torque_per_motor_Nm - ...
                   (C.drive_torque_per_motor_Nm + C.regen_torque_per_motor_Nm)));

% ---- Time base ---------------------------------------------------------
% The `time` column in this file is written with only 3 significant figures
% once t >= 1000 s ("1.00E+03"), so diff(time) collapses to 0 over the last
% ~9,980 rows and spikes to 10 s at the few places it does tick. run_sim_opt
% treats tbrake <= 0 by holding the temperature, so using diff(time) would
% silently discard the last 23% of the run. The dt_s column is the leading
% step (dt_s(i) = time(i+1) - time(i)) and agrees with diff(time) to 1e-4 s
% everywhere the time column is not saturated, so it is used instead.
dt_diff = [NaN; diff(C.time)];
dt_lead = C.dt_s;
okRegion = [false; C.time(2:end) < 1000];
V.dt_s_vs_diff_maxerr = max(abs(dt_lead(1:end-1) - dt_diff(2:end)) .* okRegion(2:end));
V.n_nonpositive_difftime = sum(dt_diff(2:end) <= 0);
V.n_difftime_gt_100ms    = sum(dt_diff(2:end) > 0.1);

if P.use_dt_s_column
    tstep = [NaN; C.dt_s(1:end-1)];          % TRAILING step for row i
    t     = [0; cumsum(C.dt_s(1:end-1))];    % reconstructed monotone clock
else
    tstep = dt_diff;
    t     = C.time;
end
V.span_used = t(end) - t(1);

% ---- Speed cleaning, exactly as compute_derived_quantities() does ------
velx = C.speed_mps;
velx(abs(velx) > P.velx_threshold) = 0;
velx(isnan(velx)) = 0;
velx(velx < 0)    = 0;
V.n_speed_cleaned = sum(velx ~= C.speed_mps);

% ---- Lap detection -----------------------------------------------------
% Each lap starts with two artifact rows (speed jumps down then back up).
% The DOWN jump is exactly the original's DS < -2.5 condition.
DS = [0; diff(velx)];
lapIdx = find(DS < -2.5);
V.lap_boundary_idx = lapIdx(:)';
V.lap_count        = numel(lapIdx) + 1;
V.lap_period_s     = median(diff(t(lapIdx)));
artifact = false(n,1);
artifact([1 2]) = true;                        % first lap's own two rows
for k = 1:numel(lapIdx)
    artifact(lapIdx(k)) = true;
    if lapIdx(k)+1 <= n, artifact(lapIdx(k)+1) = true; end
end
if ~P.exclude_lap_artifacts, artifact(:) = false; end
lapStart = false(n,1); lapStart(1) = true; lapStart(lapIdx) = true;

% ---- Wheel-radius consistency check (diagnostic only) ------------------
omega_wheel = C.motor_rpm / P.gear_ratio * 2*pi/60;
r_implied   = velx ./ max(omega_wheel, eps);
V.r_implied_mean = mean(r_implied);
V.r_implied_vs_WheelR_pct = 100*(V.r_implied_mean/P.WheelR - 1);

%% ================== PER-STEP ENERGY CHAIN (rotor-independent) ==========
% Nothing here depends on RotorMass or RotorArea, so it is computed once.
omega_all  = repmat(omega_wheel, 1, 4);          % [fl fr rl rr], identical
omega_motor = C.motor_rpm * 2*pi/60;

speed_mph = velx * 2.23694;
if P.drs_state == 1
    F_aero = P.aero_open_a*speed_mph.^2 + P.aero_open_b*speed_mph + P.aero_open_c;
else
    F_aero = P.aero_closed_a*speed_mph.^2 + P.aero_closed_b*speed_mph + P.aero_closed_c;
end
F_aero = max(F_aero, 0);

dstep = zeros(n,1);
dstep(2:end) = (velx(1:end-1) + velx(2:end))/2 .* tstep(2:end);
Edrag = zeros(n,1);
Edrag(2:end) = F_aero(2:end) .* dstep(2:end);

% Vehicle + per-corner rotational KE loss
E1 = zeros(n,1);  E1(2:end) = 0.5*P.VehicleMass*(velx(1:end-1).^2 - velx(2:end).^2);
E2 = zeros(n,1);
for c = 1:4
    E2(2:end) = E2(2:end) + 0.5*P.I_corner*(omega_all(1:end-1,c).^2 - omega_all(2:end,c).^2);
end
Energy   = E1 + E2;
AeroFrac = min(Edrag ./ max(Energy,1), 1);

% Regen: per-motor power, gated on deceleration (derived from the speed
% change, since the CSV has no accelx channel).
decel = DS < 0;
P_motor_regen = min(C.regen_torque_per_motor_Nm .* omega_motor, 0) .* double(decel);
P_veh_regen   = 4 * P_motor_regen;
ff = min(max(P.regen_front_frac_a + P.regen_front_frac_b*velx, 0), 1);
axleRegen.front = ff .* P_veh_regen;
axleRegen.rear  = (1-ff) .* P_veh_regen;

% Braking-applied flag: PEDAL FORCE ONLY. This is the sole use of
% friction_brake_force_N anywhere in this script - it never enters an
% energy term.
DSskip = DS; DSskip(DSskip < -2.5) = 0;
brakeApplied = false(n,1);
brakeApplied(2:end) = (DSskip(2:end) < 0) & ...
                      (C.friction_brake_force_N(2:end) > P.min_pedal_force_N);
brakeApplied = brakeApplied & ~artifact;

Tbias.front = repmat(P.Tbias_front,     n, 1);
Tbias.rear  = repmat(1 - P.Tbias_front, n, 1);

sides = {'front','rear'};
for s = 1:2
    sd = sides{s};
    Ecorner = Energy .* Tbias.(sd) * 0.5 .* (1 - AeroFrac);
    Eregen  = abs(axleRegen.(sd)) .* [0; tstep(2:end)] * 0.5;
    CYC.(sd).Ecorner = Ecorner;
    CYC.(sd).Eregen  = Eregen;
    CYC.(sd).Efric   = max(Ecorner - Eregen, 0);
    CYC.(sd).clamped = brakeApplied & (Ecorner - Eregen <= 0);
end
CYC.tstep   = [0; tstep(2:end)];
CYC.velx    = velx;
CYC.active  = brakeApplied;
CYC.t       = t;
CYC.lapStart= lapStart;
CYC.Energy  = Energy;
CYC.Edrag   = Edrag;

%% ================== DIAGNOSTICS ==================
act = CYC.active;
D.active_steps    = sum(act);
D.KEloss_active_kJ = sum(Energy(act))/1e3;
D.KEloss_alldecel_kJ = sum(Energy(DS<0))/1e3;
D.Edrag_active_kJ = sum(Edrag(act))/1e3;
D.regen_total_kJ  = sum(abs(P_veh_regen).*CYC.tstep)/1e3;
for s = 1:2
    sd = sides{s};
    D.(['Efric_' sd '_kJ'])   = sum(CYC.(sd).Efric(act))/1e3;
    D.(['Eregen_' sd '_kJ'])  = sum(CYC.(sd).Eregen(act))/1e3;
    D.(['clamp_' sd '_pct'])  = 100*sum(CYC.(sd).clamped)/max(sum(act),1);
end
D.friction_share_of_KE_pct = 100*(D.Efric_front_kJ + D.Efric_rear_kJ)*2 / max(D.KEloss_active_kJ,1);
D.aero_share_of_KE_pct     = 100*D.Edrag_active_kJ / max(D.KEloss_active_kJ,1);

lapEdges = [find(lapStart); n+1];
D.perLap_Efric_front_kJ = zeros(numel(lapEdges)-1,1);
D.perLap_Efric_rear_kJ  = zeros(numel(lapEdges)-1,1);
for k = 1:numel(lapEdges)-1
    rows = lapEdges(k):lapEdges(k+1)-1;
    D.perLap_Efric_front_kJ(k) = sum(CYC.front.Efric(rows) .* act(rows))/1e3;
    D.perLap_Efric_rear_kJ(k)  = sum(CYC.rear.Efric(rows)  .* act(rows))/1e3;
end
end


function print_diagnostics(V, D, P, cfg, sides, isNewFormat) %#ok<INUSL>
if isNewFormat
    fprintf('\n================ DRIVE-CYCLE VALIDATION (%s) ================\n', cfg.cycleFile);
    fprintf('rows                       : %d\n', V.rows);
    fprintf('missing values             : %d\n', V.nan_count);
    fprintf('time span                  : %.2f s   dt = %.4f s (uniform: %s)\n', ...
        V.span, V.dt_median, ternary(V.dt_uniform,'yes','no'));
    fprintf('speed values cleaned       : %d\n', V.n_speed_cleaned);
    fprintf('laps (Lap column)          : %d   (%d rows/lap, equal: %s)\n', ...
        V.lap_count, V.lap_rows_first, ternary(V.lap_rows_equal,'yes','no'));
    fprintf('implied wheel radius       : %.4f m vs WheelR %.4f m (%+.1f%%)\n', ...
        V.r_implied_mean, P.WheelR, V.r_implied_vs_WheelR_pct);
    fprintf('braking-active steps       : %d of %d  (RequiredBrakingMagnitude_N > 0)\n', ...
        V.active_steps, V.rows);
    fprintf('identity FB+AR+RB ~= 0     : max |err| = %.2e N\n', V.identity_FB_AR_RB_maxerr);
    fprintf('identity LF ~= -RB-AD (braking rows): max |err| = %.2e N\n', ...
        V.identity_LF_RB_AD_maxerr_braking);
    if V.n_active_with_AccelX_nonneg > 0
        warning('RotorSizing:UnexpectedAccel', ...
            '%d braking-active rows have AccelX_g >= 0 (expected negative during braking).', ...
            V.n_active_with_AccelX_nonneg);
    end
    fprintf('\n================ ENERGY DIAGNOSTICS ================\n');
    fprintf('required braking energy    : %.2f kJ\n', D.required_braking_total_kJ);
    fprintf('regen applied              : %.2f kJ   (available: %.2f kJ)\n', ...
        D.regen_applied_total_kJ, D.regen_available_total_kJ);
    fprintf('friction share of required : %.1f %%\n', D.friction_share_of_required_pct);
    for s = 1:2
        sd = sides{s};
        fprintf('%-5s rotor: friction %.2f kJ (direct measurement, no clamp)\n', ...
            sd, D.(['Efric_' sd '_kJ']));
    end
else
    fprintf('\n================ DRIVE-CYCLE VALIDATION (%s) ================\n', cfg.cycleFile);
    fprintf('rows                       : %d\n', V.rows);
    fprintf('missing values             : %d\n', V.nan_count);
    fprintf('time column span           : %.1f s (as written)\n', V.span_raw);
    fprintf('time base used             : %.3f s (from %s)\n', V.span_used, ...
        ternary(P.use_dt_s_column,'dt_s','diff(time)'));
    fprintf('diff(time) <= 0 rows       : %d   (>100 ms rows: %d)\n', ...
        V.n_nonpositive_difftime, V.n_difftime_gt_100ms);
    fprintf('dt_s vs diff(time), t<1000 : max |err| = %.2e s\n', V.dt_s_vs_diff_maxerr);
    fprintf('total = 4 x per-motor      : max |err| = %.2e Nm\n', V.x4_maxerr);
    fprintf('laps detected              : %d   period %.4f s\n', V.lap_count, V.lap_period_s);
    fprintf('lap boundary rows          : %s ...\n', mat2str(V.lap_boundary_idx(1:min(5,end))));
    fprintf('implied wheel radius       : %.4f m vs WheelR %.4f m (%+.1f%%)\n', ...
        V.r_implied_mean, P.WheelR, V.r_implied_vs_WheelR_pct);
    fprintf('\n================ ENERGY DIAGNOSTICS ================\n');
    fprintf('braking-active steps       : %d of %d\n', D.active_steps, V.rows);
    fprintf('KE loss on active steps    : %.1f kJ\n', D.KEloss_active_kJ);
    fprintf('  aero share               : %.1f %%\n', D.aero_share_of_KE_pct);
    fprintf('  friction share (4 rotors): %.1f %%\n', D.friction_share_of_KE_pct);
    fprintf('vehicle regen energy       : %.1f kJ\n', D.regen_total_kJ);
    for s = 1:2
        sd = sides{s};
        fprintf('%-5s rotor: friction %.1f kJ | regen %.1f kJ | per lap %.2f kJ | max(...,0) clamp on %.1f%% of braking steps\n', ...
            sd, D.(['Efric_' sd '_kJ']), D.(['Eregen_' sd '_kJ']), ...
            mean(D.(['perLap_Efric_' sd '_kJ'])), D.(['clamp_' sd '_pct']));
        if D.(['clamp_' sd '_pct']) > 5
            warning('RotorSizing:ClampActive', ...
                ['*** %s: the max(Energy_corner - regen, 0) clamp is active on %.1f%% of braking ' ...
                 'steps. Regen is being credited more energy than that corner has available, which ' ...
                 'means the regen front/rear split (regen_front_frac_*) and/or Tbias_front are ' ...
                 'inconsistent with the data. The %s heat figure is a LOWER BOUND until those are fixed.'], ...
                 upper(sd), D.(['clamp_' sd '_pct']), sd);
        end
    end
end
end


function P = load_rotor_params(xlsxPath)
%LOAD_ROTOR_PARAMS Read the shared "Parameters" sheet BY NAME into a struct.
% Mirrors loadVehicleParams.m's own import loop exactly (same read call, same
% per-row conversion, same isvarname error-on-bad-name behavior, same
% fail-loudly missing-field check) so the two loaders behave identically
% against the one spreadsheet both scripts now read.
if ~isfile(xlsxPath)
    error('RotorSizing:NoParams','Parameter file not found: %s', xlsxPath);
end
T = readtable(xlsxPath, 'Sheet', 'Parameters', 'TextType', 'string');
requiredCols = {'Name', 'Value'};
missingCols = requiredCols(~ismember(requiredCols, T.Properties.VariableNames));
if ~isempty(missingCols)
    error('RotorSizing:BadSchema', ...
        'Sheet "Parameters" is missing required column(s): %s', strjoin(missingCols, ', '));
end

P = struct();
for i = 1:height(T)
    rawName = T.Name(i);
    if ismissing(rawName)          % blank Name cell -> skip, don't error on char()
        continue
    end
    name = strtrim(char(rawName));
    if isempty(name)
        continue
    end
    if ~isvarname(name)
        error('RotorSizing:BadName', ...
            'Row %d: "%s" is not a valid MATLAB field name.', i, name);
    end
    value = T.Value(i);
    if iscell(value); value = value{1}; end
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        error('RotorSizing:BadValue', ...
            'Parameter "%s" (row %d) must be a finite numeric scalar.', name, i);
    end
    if isfield(P, name)
        error('RotorSizing:DuplicateName', ...
            'Parameter "%s" is defined more than once.', name);
    end
    P.(name) = value;
end

% ---- Presence check: fail loudly on a renamed/deleted row rather than
%      silently propagating an undefined field into the physics ----
expected = {'VehicleMass','I_corner','WheelR','gear_ratio','rho_rotor','velx_threshold', ...
    'BrakeFrac','CalibrationFactor','aero_open_a','aero_open_b','aero_open_c', ...
    'aero_closed_a','aero_closed_b','aero_closed_c','drs_state', ...
    'hw_x1_front','hw_b1_front','hw_x1_rear','hw_b1_rear','padfrac_x2','padfrac_b2', ...
    'Tbias_front','regen_front_frac_a','regen_front_frac_b','min_pedal_force_N', ...
    'T_target_F','TambC_fallback','T_init_F', ...
    'Do_front_mm','Di_front_mm','t_min_front_mm','t_max_front_mm','s_min_front_mm', ...
    'phi_max_front','cut_type_front', ...
    'Do_rear_mm','Di_rear_mm','t_min_rear_mm','t_max_rear_mm','s_min_rear_mm', ...
    'phi_max_rear','cut_type_rear','include_edge_area','edge_margin_mm', ...
    'enforce_pattern_packing','use_dt_s_column', ...
    'exclude_lap_artifacts','n_mass_grid','A_bracket_lo','A_bracket_hi', ...
    'area_bisect_iters','euler_warn_threshold'};
missingParams = expected(~isfield(P, expected));
if ~isempty(missingParams)
    error('RotorSizing:MissingParams', ...
        'Spreadsheet is missing required parameter(s): %s', strjoin(missingParams, ', '));
end
end


function g = make_geom(Do, Di, tmin, tmax, smin, phimax, cutType, includeEdges, rho, ...
                       edgeMargin, enforcePacking)
g.Do = Do; g.Di = Di; g.tmin = tmin; g.tmax = tmax; g.smin = smin;
g.rho = rho; g.edge = edgeMargin;
g.cutC = 4;  g.cutKind = 'drilled hole';
if cutType == 2, g.cutC = 2; g.cutKind = 'radial slot'; end
g.Aring  = pi/4*(Do^2 - Di^2);
g.Pedge  = 0; if includeEdges, g.Pedge = pi*(Do + Di); end
g.phimax_spec = phimax;
g.phimax_pack = phi_packing_max(g);
g.phimax = phimax;
g.packingLimited = false;
if enforcePacking && g.phimax_pack < phimax
    g.phimax = g.phimax_pack;
    g.packingLimited = true;
end
end


function phi = phi_packing_max(g)
%PHI_PACKING_MAX Largest removed-area fraction a real single-row cut pattern
% can deliver: features >= s_min, ligaments >= s_min, and a radial margin of
% edgeMargin left at both the OD and the ID. Scanned over feature size.
w  = (g.Do - g.Di)/2;
Dm = (g.Do + g.Di)/2;
featMax = w - 2*g.edge;
if featMax < g.smin, phi = 0; return; end
best = 0;
for d = linspace(g.smin, featMax, 2000)
    N = floor(pi*Dm/(d + g.smin));
    if N < 1, continue; end
    if g.cutC == 4
        Ar = N*pi*d^2/4;
    else
        Ar = N*d*(w - 2*g.edge);
    end
    best = max(best, Ar/g.Aring);
end
phi = min(best, 1);
end

function m = geo_mass(g, t, Ar),  m = g.rho * t * (g.Aring - Ar);  end
function A = geo_area(g, t, Ar, L), A = 2*(g.Aring - Ar) + L*t + g.Pedge*t; end
function L = Lmin_of(Ar), L = 2*sqrt(pi*max(Ar,0)); end
function L = Lmax_of(g, Ar), L = g.cutC*max(Ar,0)/g.smin; end

function [lo, hi, ok] = t_range_for_mass(g, m)
Vv = m/g.rho;
lo = max(g.tmin, Vv/g.Aring);
hi = g.tmax;
if g.phimax < 1, hi = min(g.tmax, Vv/((1-g.phimax)*g.Aring)); end
ok = lo <= hi*(1+1e-9) + 1e-12;
hi = max(hi, lo);
end

function [A, t, Ar, L] = geo_area_max(g, m)
%GEO_AREA_MAX Maximise area at fixed mass over (t, Ar, L<=Lmax).
% With V = m/rho fixed, Ar = Aring - V/t, so
%   A(t) = 2V/t + t*(cutC*Aring/smin + Pedge) - cutC*V/smin,
% which is convex in t: the maximum is at an endpoint of the feasible t
% range. A 41-point grid backstop is evaluated as well.
[lo, hi, ok] = t_range_for_mass(g, m);
if ~ok, A = NaN; t = NaN; Ar = NaN; L = NaN; return; end
cand = unique([lo, hi, linspace(lo, hi, 41)]);
A = -inf; t = lo; Ar = 0; L = 0;
for k = 1:numel(cand)
    tk = cand(k);
    Ark = min(max(g.Aring - (m/g.rho)/tk, 0), g.phimax*g.Aring);
    Lk  = Lmax_of(g, Ark);
    Ak  = geo_area(g, tk, Ark, Lk);
    if Ak > A, A = Ak; t = tk; Ar = Ark; L = Lk; end
end
end

function A = geo_area_min(g, m)
%GEO_AREA_MIN Minimum area at fixed mass (single cut, L = Lmin).
[lo, hi, ok] = t_range_for_mass(g, m);
if ~ok, A = NaN; return; end
A = inf;
for tk = linspace(lo, hi, 401)
    Ark = min(max(g.Aring - (m/g.rho)/tk, 0), g.phimax*g.Aring);
    A = min(A, geo_area(g, tk, Ark, Lmin_of(Ark)));
end
end

function pr = existence_proof(g, m, A_req)
%EXISTENCE_PROOF One concrete, buildable (t, Ar, L) at mass m with A >= A_req.
% Uses a SINGLE ROW of equal cuts on the mean pitch circle, sized so the
% feature and every ligament respect s_min and the radial edge margins.
% This is a proof of feasibility only - it is NOT the recommended design.
[lo, hi, ok] = t_range_for_mass(g, m);
if ~ok, error('RotorSizing:NoGeom','No geometry reproduces m = %.4f kg.', m); end
w  = (g.Do - g.Di)/2;                 % annulus radial width
Dm = (g.Do + g.Di)/2;                 % mean pitch-circle diameter
best = [];
tList = unique([lo, hi, linspace(lo, hi, 41)]);
for q = 1:numel(tList)
    t = tList(q);
    Ar_need = g.Aring - (m/g.rho)/t;
    if Ar_need < -1e-12 || Ar_need > g.phimax*g.Aring + 1e-12, continue; end
    Ar_need = min(max(Ar_need, 0), g.phimax*g.Aring);
    cand = [];
    if Ar_need <= 1e-12
        cand = struct('t',t,'Ar',0,'L',0,'N',0,'feat',0);
    elseif g.cutC == 4
        % N equal circular holes, one row on the mean pitch circle.
        % Radial edge margin >= s_min on both sides; ligament >= s_min.
        dmax = max(w - 2*g.edge, g.smin);
        for N = 1:5000
            d = 2*sqrt(Ar_need/(pi*N));             % hits Ar exactly
            if d > dmax, continue; end
            if d < g.smin, break; end
            if N*(d + g.smin) > pi*Dm, continue; end
            c = struct('t',t,'Ar',Ar_need,'L',N*pi*d,'N',N,'feat',d);
            if isempty(cand) || c.L > cand.L, cand = c; end
        end
    else
        % N radial slots of width d and radial length Lr, one row.
        Lr = max(w - 2*g.edge, g.smin);
        for N = 1:5000
            d = Ar_need/(N*Lr);                     % hits Ar exactly
            if d < g.smin, break; end
            if N*(d + g.smin) > pi*Dm, continue; end
            c = struct('t',t,'Ar',Ar_need,'L',N*2*(d+Lr),'N',N,'feat',d);
            if isempty(cand) || c.L > cand.L, cand = c; end
        end
    end
    if isempty(cand), continue; end
    cand.A = geo_area(g, cand.t, cand.Ar, cand.L);
    if cand.A >= A_req - 1e-12 && (isempty(best) || cand.A > best.A), best = cand; end
end
if isempty(best)
    error('RotorSizing:NoProof', ...
        'No single-row cut pattern at m = %.4f kg reaches A_req = %.5f m^2.', m, A_req);
end
pr = best;
pr.Dm = Dm; pr.kind = g.cutKind;
pr.m_check = geo_mass(g, pr.t, pr.Ar);
pr.A_check = geo_area(g, pr.t, pr.Ar, pr.L);
end

function A = area_frontier(peakFun, m, T_target_F, lo, hi, iters, hscale)
%AREA_FRONTIER Smallest area with peak temperature <= target, by geometric
% bisection (peak T is monotone decreasing in A).
if peakFun(m, hi, hscale) > T_target_F, A = inf; return; end
if peakFun(m, lo, hscale) <= T_target_F, A = 0; return; end
for k = 1:iters
    mid = sqrt(lo*hi);
    if peakFun(m, mid, hscale) > T_target_F, lo = mid; else, hi = mid; end
end
A = hi;
end


function TK = run_sim_rotor(CYC, sd, RotorMass, RotorArea, hwCoef, padp, P, TambK, TinitK, hscale)
%RUN_SIM_ROTOR Thermal march, identical to run_sim_opt() in
% BrakeCoeffOptimizer.m. Deviations are confined to the INPUTS:
%   - the braking flag comes from pedal force, not brake pressure;
%   - regen comes from motor torque x motor rpm;
%   - omega comes from motor_rpm / gear_ratio;
%   - the per-step energy chain up to friction_energy is precomputed in CYC
%     (same formulae, same order, same max(...,0) clamp).
% The branch logic, the order of operations (heat added first, convection
% then computed from the UPDATED temperature) and the temperature-dependent
% SpecHeat are unchanged.
tstep = CYC.tstep;  velx = CYC.velx;  act = CYC.active;
Efric = CYC.(sd).Efric;
nSteps = numel(tstep);
TK = zeros(nSteps,1);  TK(1) = TinitK;
BFCF = P.BrakeFrac * P.CalibrationFactor;
for i = 2:nSteps
    prevTemp = TK(i-1);
    tbrake   = tstep(i);
    if tbrake <= 0, TK(i) = prevTemp; continue; end

    h_w      = (hwCoef(1)*velx(i) + hwCoef(2)) * hscale;
    Rrotor   = 1/(h_w*RotorArea);
    PadFrac  = max(min(padp(1)*prevTemp + padp(2), 1), 0);
    SpecHeat = (0.0005*prevTemp + 0.2813)*1000;

    if act(i)
        CorrectedEnergy = Efric(i) * (1 - PadFrac) * BFCF;
        TK(i) = CorrectedEnergy/(RotorMass*SpecHeat) + prevTemp;
        qout  = (TK(i) - TambK)/Rrotor;
        TK(i) = TK(i) - (qout*tbrake)/(RotorMass*SpecHeat);
    else
        qout  = (prevTemp - TambK)/Rrotor;
        TK(i) = prevTemp - (qout*tbrake)/(RotorMass*SpecHeat);
    end
end
end

function F = K2F(K), F = (K - 273.15)*(9/5) + 32; end
function out = ternary(c, a, b), if c, out = a; else, out = b; end, end
