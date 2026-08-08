%% ================================================================
%  Brake Model Coefficient Optimizer
% ================================================================
%  Fits, across MULTIPLE datasets simultaneously:
%    - Rotor cooling coefficient:  h_w(v) = x1*v + b1
%    - Pad energy fraction:        PadFrac = f(T_rotor, P_applied)
%
%  Several candidate forms for PadFrac(T,P) are tried and compared
%  (see CANDIDATE PAD-FRACTION MODELS below) so you can pick the one
%  that generalizes best, not just the one with the lowest raw error.
%
%  Requires: Optimization Toolbox (lsqnonlin)
%
%  Datasets are selected via a multi-select file picker. Each dataset
%  can have its own time crop and its own column layout (21/22/27 col
%  formats are auto-detected from the file itself).
% ================================================================

clc; clear; close all
cfg.SkipTimeCropPrompt = true;

%% ================== USER INPUTS ==================
% Vehicle / rotor parameters (same as BrakeDataAnalysis.m - adjust as needed)
VehicleMass     = 259;      % kg
RotorMass_front = 0.5107;   % kg
RotorMass_rear  = 0.3343;   % kg
RotorArea_front = 0.038;    % m^2
RotorArea_rear  = 0.0226;   % m^2
I               = 0.30754;  % rotational inertia kg*m^2
WheelR          = 0.213;    % meters
gear_ratio      = 12.97;    % motor to wheel gear ratio
TambC           = 25;       % degC - fallback only; each dataset now derives its own ambient (see loading loop)

% Brake geometry (needed to reconstruct Tbias_brake from raw pressures)
front_piston_count = 6;
front_piston_dia   = 0.0157226;  % meters
front_rotor_dia    = 0.18542;    % meters
rear_piston_count  = 4;
rear_piston_dia    = 0.0141732;  % meters
rear_rotor_dia     = 0.18796;    % meters

% Mu vs temperature lookup table
mu_temp_table = [100, 200, 300, 400, 500, 600, 700, 800, 900, 950, 1100, 1200];  % degF
mu_table      = [0.45, 0.46, 0.49, 0.53, 0.55, 0.56, 0.565, 0.565, 0.57, 0.57, 0.535, 0.535];

% Aero drag curve coefficients (drag force in N vs speed in mph)
aero_open_a   =  0.0793024;   aero_open_b   =  1.46483;    aero_open_c   = -37.23993;
aero_closed_a =  0.136247;    aero_closed_b =  2.53449;    aero_closed_c = -67.44975;

% Velocity cleaning threshold
velx_threshold = 34;  % m/s

% Seed values (team's current fit) - used as optimizer starting point
x1_seed = 2.0006;    % h_w slope
b1_seed = 40.7137;   % h_w intercept
x2_seed = 0.001005;  % PadFrac vs temperature (Kelvin) slope
b2_seed = -0.5385;   % PadFrac intercept

%% ================== SELECT DATASETS ==================
[files, path] = uigetfile('*.txt', 'Select dataset file(s)', 'MultiSelect', 'on');
if isequal(files, 0)
    error('No files selected.');
end
if ischar(files)
    files = {files};   % single file selected -> wrap in cell for uniform handling
end
nFiles = numel(files);

aeroParams = struct('front_piston_count', front_piston_count, 'front_piston_dia', front_piston_dia, ...
    'rear_piston_count', rear_piston_count, 'rear_piston_dia', rear_piston_dia, ...
    'front_rotor_dia', front_rotor_dia, 'rear_rotor_dia', rear_rotor_dia, ...
    'mu_temp_table', mu_temp_table, 'mu_table', mu_table, ...
    'aero_open_a', aero_open_a, 'aero_open_b', aero_open_b, 'aero_open_c', aero_open_c, ...
    'aero_closed_a', aero_closed_a, 'aero_closed_b', aero_closed_b, 'aero_closed_c', aero_closed_c);
gearParams = struct('gear_ratio', gear_ratio, 'velx_threshold', velx_threshold);

datasets = struct('name', {}, 't', {}, 'velx', {}, 'frontpressure', {}, 'rearpressure', {}, ...
    'fr_temp_F', {}, 'rr_temp_F', {}, 'Tbias_brake', {}, 'F_aero', {}, 'Edrag', {}, ...
    'total_regen_power', {}, 'fl_omega_wheel', {}, 'fr_omega_wheel', {}, ...
    'rl_omega_wheel', {}, 'rr_omega_wheel', {}, 'TambK', {});

for k = 1:nFiles
    fname = fullfile(path, files{k});
    fprintf('\nLoading dataset %d/%d: %s\n', k, nFiles, files{k});

    raw = readmatrix(fname, 'FileType', 'text', 'Delimiter', '\t', ...
        'NumHeaderLines', 6, 'TreatAsMissing', {'', 'NaN'});

    ncols = size(raw, 2);
    switch ncols
        case 27, fmt = 'drs_imu';
        case 22, fmt = 'drs';
        case 21, fmt = 'no_drs';
        otherwise
            error('Unrecognized column count (%d) in %s. Expected 21, 22, or 27.', ncols, files{k});
    end
    fprintf('  Detected format: %s (%d columns)\n', fmt, ncols);

    % Per-dataset time crop
    t_full = raw(:, 1);
    
    if cfg.SkipTimeCropPrompt
        t_start_k = t_full(1);
        t_end_k   = t_full(end);
    else
        prompt = {sprintf('Start time (s)  [data range %.1f - %.1f]:', t_full(1), t_full(end)), 'End time (s):'};
        defaultAns = {num2str(t_full(1)), num2str(t_full(end))};
        answer = inputdlg(prompt, sprintf('Time crop: %s', files{k}), 1, defaultAns);
        if isempty(answer)
            t_start_k = t_full(1); t_end_k = t_full(end);
        else
            t_start_k = str2double(answer{1});
            t_end_k   = str2double(answer{2});
        end
    end

    raw = raw(t_full >= t_start_k & t_full <= t_end_k, :);
    if size(raw, 1) < 3
        warning('  Dataset %s has < 3 samples after crop - skipping.', files{k});
        continue
    end

    parsed  = parse_dataset_columns(raw, fmt);
    derived = compute_derived_quantities(parsed, gearParams, aeroParams);
    Edrag   = compute_edrag(parsed.t, derived.velx, derived.F_aero);

    % Per-dataset ambient temperature: average of this dataset's first
    % front and rear rotor-temperature readings (degF -> degC -> K).
    % Falls back to the global TambC if a first reading is missing.
    firstFrontF = derived.fr_temp_F(1);
    firstRearF  = derived.rr_temp_F(1);
    if isfinite(firstFrontF) && isfinite(firstRearF)
        TambC_dataset = ((firstFrontF + firstRearF) / 2 - 32) * (5/9);
    elseif isfinite(firstFrontF)
        TambC_dataset = (firstFrontF - 32) * (5/9);
        warning('  Dataset %s: first rear temperature reading missing - ambient taken from front only.', files{k});
    elseif isfinite(firstRearF)
        TambC_dataset = (firstRearF - 32) * (5/9);
        warning('  Dataset %s: first front temperature reading missing - ambient taken from rear only.', files{k});
    else
        TambC_dataset = TambC;
        warning('  Dataset %s: first front/rear temperature readings both missing - falling back to global TambC = %.1f C.', files{k}, TambC);
    end
    TambK_dataset = TambC_dataset + 273.15;
    fprintf('  Dataset ambient temperature: %.1f C (%.1f F)\n', TambC_dataset, TambC_dataset*9/5+32);

    idx = numel(datasets) + 1;
    datasets(idx).name              = files{k};
    datasets(idx).t                 = parsed.t;
    datasets(idx).velx              = derived.velx;
    datasets(idx).frontpressure     = derived.frontpressure;
    datasets(idx).rearpressure      = derived.rearpressure;
    datasets(idx).fr_temp_F         = derived.fr_temp_F;
    datasets(idx).rr_temp_F         = derived.rr_temp_F;
    datasets(idx).Tbias_brake       = derived.Tbias_brake;
    datasets(idx).F_aero            = derived.F_aero;
    datasets(idx).Edrag             = Edrag;
    datasets(idx).total_regen_power = derived.total_regen_power;
    datasets(idx).fl_omega_wheel    = derived.fl_omega_wheel;
    datasets(idx).fr_omega_wheel    = derived.fr_omega_wheel;
    datasets(idx).rl_omega_wheel    = derived.rl_omega_wheel;
    datasets(idx).rr_omega_wheel    = derived.rr_omega_wheel;
    datasets(idx).TambK             = TambK_dataset;

    fprintf('  %d samples, %.1f s duration\n', numel(parsed.t), parsed.t(end) - parsed.t(1));
end
nFiles = numel(datasets);
if nFiles == 0
    error('No usable datasets loaded.');
end

%% ================== DATA-DRIVEN COEFFICIENT BOUNDS ==================
% Bounds are derived from the actual observed speed/temperature/pressure
% ranges across all loaded datasets, so they scale to physical reality
% instead of being arbitrary fixed numbers.
allT_K = []; allP = [];
for k = 1:nFiles
    allT_K = [allT_K; (datasets(k).fr_temp_F - 32)*(5/9) + 273.15; ...
                       (datasets(k).rr_temp_F - 32)*(5/9) + 273.15]; %#ok<AGROW>
    allP   = [allP; datasets(k).frontpressure; datasets(k).rearpressure]; %#ok<AGROW>
end
Tmin_K = min(allT_K); Tmax_K = max(allT_K);
dT = max(Tmax_K - Tmin_K, 1);
Pmax = max(allP);
dP = max(Pmax, 1);

% h_w = x1*v + b1  [W/m^2K]. Must stay strictly positive (convection
% coefficient can't be zero/negative). Bounded around the team's
% existing fit as an informed prior - widen if you don't trust that prior.
x1_lb = 0;      x1_ub = 3 * x1_seed;
b1_lb = 5;      b1_ub = 3 * b1_seed;   % b1_lb=5 keeps h_w>0 at v=0

% PadFrac coefficient bounds: sized so each term, swept across the FULL
% observed T/P range, can move PadFrac by at most ~1.5 in magnitude.
% This is intentionally loose (PadFrac is hard-clamped to [0,1] at
% simulation runtime regardless) - its only job is to keep the search
% well-conditioned instead of letting coefficients blow up.
x2_bound = 1.5 / dT;          % 1/K
x3_bound = 1.5 / dP;          % 1/psi
x4_bound = 1.5 / (dT * dP);   % 1/(K*psi)
b2_bound = 1.5;

fprintf('\nData ranges used for bounds: T = [%.1f, %.1f] K, P_max = %.1f psi\n', Tmin_K, Tmax_K, Pmax);
Tmid_K = (Tmin_K + Tmax_K) / 2;
Pmid   = Pmax / 2;

%% ================== PHYSICALLY-GROUNDED PADFRAC ANCHOR ==================
% Classical heat-partition theory (Charron's relation / Newcomb) gives the
% steady-state fraction of frictional energy absorbed by each body from
% their thermal effusivity xi = sqrt(k*rho*cp), for equal/overlapping
% contact area (the pad's footprint on the rotor is the same patch for
% both bodies, so area cancels out of the ratio):
%
%   PadFrac_ideal = xi_pad / (xi_pad + xi_rotor)
%
% Unlike a pressure term in an unconstrained linear fit, this value
% cannot run away to 0 or blow up with extrapolation - it's a fixed
% material property ratio. Physically, published thermal-contact-
% conductance literature (Cooper-Mikic-Yovanovich and related work)
% shows contact conductance INCREASING (not decreasing) with clamping
% pressure as surface asperities flatten, so if anything real pressure
% dependence should pull PadFrac TOWARD this ideal value as pressure
% rises, not away from it - the opposite of collapsing to zero.

% Rotor: normalized 4130 alloy steel (MatWeb datasheet, ~100-300 C range)
k_rotor   = 42.7;    % W/m-K
rho_rotor = 7850;    % kg/m^3
cp_rotor  = 477;     % J/kg-K
xi_rotor  = sqrt(k_rotor * rho_rotor * cp_rotor);

% Pad: Porterfield R4-1 (carbon-Kevlar semi-metallic composite). Porterfield
% does not publish k/rho/cp for this compound, so these are literature
% estimates for comparable semi-metallic friction composites, NOT a
% manufacturer-certified value - treat PadFrac_ideal as a sanity-check
% anchor/starting point, not ground truth:
%   k  ~ 2-5 W/m-K   (semi-metallic composites w/ steel/Kevlar fiber content)
%   rho~ 2000-2500 kg/m^3
%   cp ~ 800-1200 J/kg-K
k_pad_lo = 2.0; k_pad_mid = 3.0; k_pad_hi = 5.0;
rho_pad_lo = 2000; rho_pad_mid = 2200; rho_pad_hi = 2500;
cp_pad_lo = 800; cp_pad_mid = 1000; cp_pad_hi = 1200;

xi_pad_lo  = sqrt(k_pad_lo  * rho_pad_lo  * cp_pad_lo);
xi_pad_mid = sqrt(k_pad_mid * rho_pad_mid * cp_pad_mid);
xi_pad_hi  = sqrt(k_pad_hi  * rho_pad_hi  * cp_pad_hi);

PadFrac_ideal_lo  = xi_pad_lo  / (xi_pad_lo  + xi_rotor);
PadFrac_ideal     = xi_pad_mid / (xi_pad_mid + xi_rotor);   % central estimate, used below
PadFrac_ideal_hi  = xi_pad_hi  / (xi_pad_hi  + xi_rotor);

fprintf('\n================ PHYSICALLY-GROUNDED PADFRAC ANCHOR ================\n');
fprintf('Rotor (4130 steel) effusivity:  %.0f J/(m^2*K*sqrt(s))\n', xi_rotor);
fprintf('Pad (R4-1, estimated) effusivity: %.0f J/(m^2*K*sqrt(s)) [range %.0f-%.0f]\n', ...
    xi_pad_mid, xi_pad_lo, xi_pad_hi);
fprintf('PadFrac_ideal = %.3f  (estimated range %.3f - %.3f)\n', ...
    PadFrac_ideal, PadFrac_ideal_lo, PadFrac_ideal_hi);
fprintf('This is a starting-point anchor, not a certified value - Porterfield does\n');
fprintf('not publish R4-1''s thermal properties. Used to seed/bound models 5 and 6 below.\n');

%% ================== CANDIDATE PAD-FRACTION MODELS ==================
% Each model computes the RAW (pre-clamp) PadFrac from rotor temp (K)
% and applied pressure (psi). Runtime clamping to [0,1] happens inside
% run_sim_opt, exactly as in the original brake_temp_sim.m.
models = struct('name', {}, 'fun', {}, 'lb', {}, 'ub', {}, 'x0', {});

models(1).name = 'Linear in T only (current baseline)';
models(1).fun  = @(T, P, p) p(1).*T + p(2);
models(1).lb   = [-x2_bound, -b2_bound];
models(1).ub   = [ x2_bound,  b2_bound];
models(1).x0   = [x2_seed, b2_seed];

models(2).name = 'Linear, independent T and P';
models(2).fun  = @(T, P, p) p(1).*T + p(2).*P + p(3);
% p(2) (direct pressure term) is bounded >= 0: thermal-contact-conductance
% literature (Cooper-Mikic-Yovanovich) shows contact conductance can only
% increase or stay flat with clamping pressure as asperities flatten, never
% decrease - a negative pressure sensitivity has no physical basis, and
% allowing it is what let this model collapse PadFrac to zero at high P.
models(2).lb   = [-x2_bound, 0, -b2_bound];
models(2).ub   = [ x2_bound, x3_bound, b2_bound];
models(2).x0   = [x2_seed, 0, b2_seed];

models(3).name = 'Linear with T*P interaction';
models(3).fun  = @(T, P, p) p(1).*T + p(2).*P + p(3).*T.*P + p(4);
% p(2) (direct P term) non-negative for the same reason as model 2. The
% interaction term p(3) is left unconstrained in sign since there isn't a
% clear physical prior for how P's effect should change with T.
models(3).lb   = [-x2_bound, 0, -x4_bound, -b2_bound];
models(3).ub   = [ x2_bound, x3_bound, x4_bound,  b2_bound];
models(3).x0   = [x2_seed, 0, 0, b2_seed];

models(4).name = 'Quadratic in T, linear in P';
quad_bound = x2_bound / dT;
models(4).fun  = @(T, P, p) p(1).*T + p(2).*T.^2 + p(3).*P + p(4);
models(4).lb   = [-x2_bound, -quad_bound, 0, -b2_bound];
models(4).ub   = [ x2_bound,  quad_bound, x3_bound,  b2_bound];
models(4).x0   = [x2_seed, 0, 0, b2_seed];

models(5).name = 'Anchored to effusivity-based ideal (small T,P correction)';
% PadFrac = PadFrac_ideal + a small correction. Unlike models 1-4, this
% model doesn't fit the partition ratio from scratch - it stays tethered
% to the physically-computed anchor and only lets the optimizer nudge it
% within a tight band, which is the "decouple from the anchor" approach
% suggested by the effusivity theory above. p(1),p(2) are the correction's
% full swing across the observed T/P range; p(3) is a small constant offset.
anchor_correction_bound = 0.15;   % max additional swing away from the anchor
models(5).fun  = @(T, P, p) PadFrac_ideal + p(1).*(T-Tmid_K)/dT + p(2).*(P-Pmid)/dP + p(3);
models(5).lb   = [-anchor_correction_bound, -anchor_correction_bound, -0.05];
models(5).ub   = [ anchor_correction_bound,  anchor_correction_bound,  0.05];
models(5).x0   = [0, 0, 0];

models(6).name = 'Saturating logistic in T,P (bounded, cannot collapse to 0)';
% PadFrac = PadMax / (1 + exp(-(slope terms))). This form is bounded in
% [0, PadMax] by construction for ANY coefficient values - it cannot be
% driven to 0 or blow up by extrapolating past the well-sampled pressure
% range the way an unbounded linear/polynomial term can. PadMax is itself
% fit, with bounds informed by the effusivity-based anchor above (given
% generous headroom since the anchor is a rough estimate).
padmax_lb = max(0.02, PadFrac_ideal_lo * 0.5);
padmax_ub = min(0.6,  PadFrac_ideal_hi * 2.5);
models(6).fun  = @(T, P, p) p(1) ./ (1 + exp(-(p(2).*(T-Tmid_K)/dT + p(3).*(P-Pmid)/dP + p(4))));
models(6).lb   = [padmax_lb, -15, -15, -10];
models(6).ub   = [padmax_ub,  15,  15,  10];
models(6).x0   = [min(max(PadFrac_ideal*1.3, padmax_lb), padmax_ub), 0, 0, 0];

%% ================== RUN OPTIMIZATION FOR EACH MODEL ==================
opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 500, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

results = struct('name', {}, 'params', {}, 'h_wF', {}, 'h_wR', {}, 'padfrac_params', {}, 'fun', {}, ...
    'rmse_F', {}, 'avg_pct_err', {}, 'sse', {}, 'nresid', {}, 'nparams', {}, 'aicc', {});

% Two fitting strategies are compared for every PadFrac model:
%   'joint'    - the original approach: h_w and PadFrac are fit together
%                against the full residual in one lsqnonlin call.
%   'twostage' - h_w is fit FIRST, using only cooldown-phase residuals
%                (steps where PadFrac plays no role at all, since no
%                energy is being added), with PadFrac held at the fixed
%                physical anchor. h_w is then FROZEN and PadFrac is fit
%                against the full residual. This breaks the confounding
%                between the two parameter sets described in the earlier
%                discussion of why the pressure terms were misbehaving:
%                in the joint fit, h_w and PadFrac can trade off against
%                each other to explain the same temperature residual,
%                which is a major source of the pressure term drifting to
%                physically implausible values.
strategies = {'joint', 'twostage'};
constPadFracFun = @(T, P, p) PadFrac_ideal; %#ok<NASGU> % placeholder for stage 1 (p unused)

idx = 0;
for m = 1:numel(models)
    for s = 1:numel(strategies)
        strategy = strategies{s};
        idx = idx + 1;
        label = sprintf('%s [%s]', models(m).name, strategy);
        fprintf('\n=== Fitting model %d/%d (%s): %s ===\n', m, numel(models), strategy, models(m).name);

        if strcmp(strategy, 'joint')
            x0 = [x1_seed, b1_seed, x1_seed, b1_seed, models(m).x0];
            lb = [x1_lb, b1_lb, x1_lb, b1_lb, models(m).lb];
            ub = [x1_ub, b1_ub, x1_ub, b1_ub, models(m).ub];

            resFun = @(p) brake_temp_residuals(p, models(m).fun, datasets, ...
                VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
                I, WheelR, 'all');

            [p_opt, ~, residual] = lsqnonlin(resFun, x0, lb, ub, opts);
        else
            % Stage 1: h_w only, cooldown-phase residual, PadFrac fixed at anchor.
            x0_hw = [x1_seed, b1_seed, x1_seed, b1_seed];
            lb_hw = [x1_lb, b1_lb, x1_lb, b1_lb];
            ub_hw = [x1_ub, b1_ub, x1_ub, b1_ub];
            resFun_stage1 = @(hw) brake_temp_residuals([hw(:); 0], constPadFracFun, datasets, ...
                VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
                I, WheelR, 'cooldown_only');
            hw_opt = lsqnonlin(resFun_stage1, x0_hw, lb_hw, ub_hw, opts);

            % Stage 2: PadFrac params only, full residual, h_w frozen at Stage-1 result.
            resFun_stage2 = @(pp) brake_temp_residuals([hw_opt(:); pp(:)], models(m).fun, datasets, ...
                VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
                I, WheelR, 'all');
            [pp_opt, ~, residual] = lsqnonlin(resFun_stage2, models(m).x0, models(m).lb, models(m).ub, opts);

            p_opt = [hw_opt(:); pp_opt(:)];
        end

        residual = residual(isfinite(residual));
        nResid   = numel(residual);
        sse      = sum(residual.^2);
        rmse     = sqrt(sse / nResid);
        kParams  = numel(p_opt);
        % AICc: lower is better; penalizes extra parameters so a model isn't
        % favored just because it has more knobs to turn.
        aicc = nResid*log(sse/nResid) + 2*kParams + (2*kParams*(kParams+1)) / max(nResid - kParams - 1, 1);

        results(idx).name           = label;
        results(idx).params         = p_opt;
        results(idx).h_wF           = p_opt(1:2);
        results(idx).h_wR           = p_opt(3:4);
        results(idx).padfrac_params = p_opt(5:end);
        results(idx).fun            = models(m).fun;
        results(idx).rmse_F         = rmse;
        results(idx).sse            = sse;
        results(idx).nresid         = nResid;
        results(idx).nparams        = kParams;
        results(idx).aicc           = aicc;

        fprintf('  h_wF:  x1 = %.4f, b1 = %.4f\n', p_opt(1), p_opt(2));
        fprintf('  h_wR:  x1 = %.4f, b1 = %.4f\n', p_opt(3), p_opt(4));
        fprintf('  PadFrac params: %s\n', mat2str(p_opt(5:end), 5));
        fprintf('  RMSE = %.2f degF | AICc = %.1f\n', rmse, aicc);
    end
end

%% ================== COMPARE MODELS ==================
fprintf('\n================ MODEL COMPARISON ================\n');
fprintf('%-58s %8s %8s %8s\n', 'Model [strategy]', '#params', 'RMSE(F)', 'AICc');
for m = 1:numel(results)
    fprintf('%-58s %8d %8.2f %8.1f\n', results(m).name, results(m).nparams, results(m).rmse_F, results(m).aicc);
end
[~, best_idx] = min([results.aicc]);
fprintf('\n>>> Recommended model (lowest AICc, balances fit vs. complexity): %s <<<\n', results(best_idx).name);
fprintf('    If you just want the lowest raw error regardless of overfit risk, compare RMSE instead.\n');

%% ================== PLOT BEST MODEL: MEASURED VS SIMULATED ==================
best = results(best_idx);
best_fun = best.fun;

for k = 1:nFiles
    ds = datasets(k);
    predF_front = run_sim_opt(ds.t, ds.velx, ds.frontpressure, ds.fr_temp_F, ds.Tbias_brake, ...
        best.h_wF(1), best.h_wF(2), best_fun, best.padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.fl_omega_wheel, ds.fr_omega_wheel, VehicleMass, RotorMass_front, RotorArea_front, I, WheelR, ds.TambK);
    predF_rear = run_sim_opt(ds.t, ds.velx, ds.rearpressure, ds.rr_temp_F, 1 - ds.Tbias_brake, ...
        best.h_wR(1), best.h_wR(2), best_fun, best.padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.rl_omega_wheel, ds.rr_omega_wheel, VehicleMass, RotorMass_rear, RotorArea_rear, I, WheelR, ds.TambK);

    figure('Name', sprintf('Best Fit - %s', ds.name));
    subplot(2,1,1);
    plot(ds.t, ds.fr_temp_F, 'k-', 'DisplayName', 'Measured Front'); hold on;
    plot(ds.t, predF_front, 'r--', 'DisplayName', 'Simulated Front');
    xlabel('Time (s)'); ylabel('Temp (degF)'); legend('Location','best'); grid on;
    title(sprintf('%s - Front Rotor', ds.name), 'Interpreter', 'none');

    subplot(2,1,2);
    plot(ds.t, ds.rr_temp_F, 'k-', 'DisplayName', 'Measured Rear'); hold on;
    plot(ds.t, predF_rear, 'b--', 'DisplayName', 'Simulated Rear');
    xlabel('Time (s)'); ylabel('Temp (degF)'); legend('Location','best'); grid on;
    title(sprintf('%s - Rear Rotor', ds.name), 'Interpreter', 'none');
end

%% ================== FINAL COEFFICIENTS ==================
fprintf('\n================ FINAL FIT (best model) ================\n');
fprintf('Date and time %s\n', datetime);
fprintf('x1f = %.6f;   %% h_wF slope\n', best.h_wF(1));
fprintf('b1f = %.6f;   %% h_wF intercept\n', best.h_wF(2));
fprintf('x1r = %.6f;   %% h_wR slope\n', best.h_wR(1));
fprintf('b1r = %.6f;   %% h_wR intercept\n', best.h_wR(2));
fprintf('padFrac function: %s\n', func2str(best_fun));
fprintf('PadFrac params (%s):\n', best.name);
for i = 1:numel(best.padfrac_params)
    fprintf('  p(%d) = %.8g\n', i, best.padfrac_params(i));
end
fprintf('\nNOTE: if the best model has more than 2 PadFrac params (i.e. is not\n');
fprintf('the plain linear-in-T baseline), you must update the PadFrac formula\n');
fprintf('inside run_sim (in brake_temp_sim.m) to match the "%s" form above -\n', best.name);
fprintf('the old hardcoded "PadFrac = prevTemp*x2 + b2" line will not use the\n');
fprintf('new pressure-dependent or logistic/anchored terms.\n');

%% ================== AGGREGATED PLOTS (Rotor temp in °F, axes start at 0) ==================
% Creates two figures total combining all parsed datasets.
% Assumes variables: datasets, best, best_fun, results, best_idx exist.

% -- prepare global ranges --
all_vel = vertcat(datasets.velx);
vmin = 0; vmax = max(all_vel);
vvec = linspace(vmin, vmax, 400);

all_TF_pts = vertcat(datasets.fr_temp_F, datasets.rr_temp_F); % measured rotor temps in °F
Tmin_F = min(all_TF_pts); Tmax_F = max(all_TF_pts);

all_P = vertcat(datasets.frontpressure, datasets.rearpressure);
Pmin = 0; Pmax = max(all_P);

% --- FIGURE 1: h_w vs vehicle velocity (all datasets combined) ---
figure('Name','h\_w vs Velocity - All Datasets');
hold on;
for k = 1:numel(datasets)
    % use same representative linear model for each dataset (over global vvec)
    hF = best.h_wF(1) .* vvec + best.h_wF(2);
    hR = best.h_wR(1) .* vvec + best.h_wR(2);
    plot(vvec, hF, 'r-', 'LineWidth', 1, 'HandleVisibility','off');
    plot(vvec, hR, 'b--', 'LineWidth', 1, 'HandleVisibility','off');
end
plot(vvec, best.h_wF(1).*vvec + best.h_wF(2), 'r-', 'LineWidth', 2, 'DisplayName','Front h\_w (all)');
plot(vvec, best.h_wR(1).*vvec + best.h_wR(2), 'b--', 'LineWidth', 2, 'DisplayName','Rear h\_w (all)');
xlabel('Vehicle Velocity (m/s)');
ylabel('h\_w (W/m^2K)');
xlim([0, vmax]);
legend('Location','best');
grid on;
title('Convective Coefficient h\_w vs Velocity (All Datasets)');

% --- FIGURE 2: PadFrac(T_{°F}, P) 3D surface with T in °F and axes starting at 0 ---
% Create TF and P grids (TF in °F), convert to K for model evaluation
Tvec_F = linspace(max(0, Tmin_F - 5), Tmax_F + 5, 160);
Pvec   = linspace(Pmin, Pmax + 2, 120);
[TFF, PP] = meshgrid(Tvec_F, Pvec);

% convert TF grid to K for model
TK_grid = (TFF - 32) * (5/9) + 273.15;
PF_raw = best_fun(TK_grid, PP, best.padfrac_params);
PF = min(max(PF_raw, 0), 1);

figure('Name','PadFrac vs Temperature(°F) and Pressure - All Datasets');
surf(TFF, PP, PF, 'EdgeColor','none', 'FaceAlpha', 0.9);
hold on;
% overlay measured sample points (TF, P)
all_TF_pts = []; all_P_pts = [];
for k = 1:numel(datasets)
    ds = datasets(k);
    TF_fr = ds.fr_temp_F(:);
    TF_rr = ds.rr_temp_F(:);
    all_TF_pts = [all_TF_pts; TF_fr; TF_rr]; %#ok<AGROW>
    all_P_pts  = [all_P_pts; ds.frontpressure(:); ds.rearpressure(:)]; %#ok<AGROW>
end
TK_pts = (all_TF_pts - 32) * (5/9) + 273.15;
PF_pts_raw = best_fun(TK_pts, all_P_pts, best.padfrac_params);
PF_pts = min(max(PF_pts_raw, 0), 1);
scatter3(all_TF_pts, all_P_pts, PF_pts, 18, PF_pts, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.85);

colorbar;
xlabel('Rotor Temp (°F)');
ylabel('Applied Pressure (psi)');
zlabel('PadFrac (clamped)');
xlim([0, max(Tvec_F)]);
ylim([0, max(Pvec)]);
title(sprintf('PadFrac(T_{°F},P) — All Datasets (model: %s)', results(best_idx).name), 'Interpreter', 'none');
view(45,25);
grid on;


%% ================================================================
%  LOCAL FUNCTIONS
% ================================================================

function parsed = parse_dataset_columns(data, fmt)
% Mirrors the "Data Parsing" section of BrakeDataAnalysis.m
parsed.t           = data(:,1);
parsed.fl_temp_adc = data(:,2);
parsed.fr_temp_adc = data(:,3);
parsed.rl_temp_adc = data(:,4);
parsed.rr_temp_adc = data(:,5);

switch fmt
    case 'drs_imu'
        parsed.drs_state         = data(:,6);
        parsed.frontpressure_adc = data(:,8);
        parsed.rearpressure_adc  = data(:,12);
        parsed.fl_Tmotor_Mn = data(:,13); parsed.fl_vmotor = data(:,14);
        parsed.fr_Tmotor_Mn = data(:,15); parsed.fr_vmotor = data(:,16);
        parsed.rl_Tmotor_Mn = data(:,17); parsed.rl_vmotor = data(:,18);
        parsed.rr_Tmotor_Mn = data(:,19); parsed.rr_vmotor = data(:,20);
        parsed.velx   = data(:,26);
        parsed.accelx = data(:,27);
    case 'drs'
        parsed.drs_state         = data(:,6);
        parsed.frontpressure_adc = data(:,8);
        parsed.rearpressure_adc  = data(:,12);
        parsed.fl_Tmotor_Mn = data(:,13); parsed.fl_vmotor = data(:,14);
        parsed.fr_Tmotor_Mn = data(:,15); parsed.fr_vmotor = data(:,16);
        parsed.rl_Tmotor_Mn = data(:,17); parsed.rl_vmotor = data(:,18);
        parsed.rr_Tmotor_Mn = data(:,19); parsed.rr_vmotor = data(:,20);
        parsed.velx   = data(:,21);
        parsed.accelx = data(:,22);
    case 'no_drs'
        parsed.drs_state         = zeros(size(data,1), 1);  % default DRS closed
        parsed.frontpressure_adc = data(:,7);
        parsed.rearpressure_adc  = data(:,11);
        parsed.fl_Tmotor_Mn = data(:,12); parsed.fl_vmotor = data(:,13);
        parsed.fr_Tmotor_Mn = data(:,14); parsed.fr_vmotor = data(:,15);
        parsed.rl_Tmotor_Mn = data(:,16); parsed.rl_vmotor = data(:,17);
        parsed.rr_Tmotor_Mn = data(:,18); parsed.rr_vmotor = data(:,19);
        parsed.velx   = data(:,20);
        parsed.accelx = data(:,21);
    otherwise
        error('Unknown dataset format: %s', fmt);
end
end


function derived = compute_derived_quantities(parsed, gearParams, aeroParams)
% Reconstructs everything brake_temp_sim needs: cleaned velocity,
% front/rear pressure (psi), front/rear rotor temp (degF), Tbias_brake,
% aero drag force, total regen power, and per-corner wheel omega.

velx = parsed.velx;
velx(abs(velx) > gearParams.velx_threshold) = 0;
velx(isnan(velx)) = 0;
velx(velx < 0)    = 0;
speed_mph = velx * 2.23694;

frontpressure = max(0.924 * parsed.frontpressure_adc - 332.64,  0);  % psi
rearpressure  = max(0.924 * parsed.rearpressure_adc  - 376.068, 0);  % psi

fr_temp_C = 0.246 * (parsed.fr_temp_adc - 406);
rr_temp_C = 0.246 * (parsed.rr_temp_adc - 406);
fr_temp_F = fr_temp_C * (9/5) + 32;
rr_temp_F = rr_temp_C * (9/5) + 32;

fl_Tmotor = parsed.fl_Tmotor_Mn / 100 * 9.8;  fr_Tmotor = parsed.fr_Tmotor_Mn / 100 * 9.8;
rl_Tmotor = parsed.rl_Tmotor_Mn / 100 * 9.8;  rr_Tmotor = parsed.rr_Tmotor_Mn / 100 * 9.8;

fl_vwheel = parsed.fl_vmotor / gearParams.gear_ratio;  fr_vwheel = parsed.fr_vmotor / gearParams.gear_ratio;
rl_vwheel = parsed.rl_vmotor / gearParams.gear_ratio;  rr_vwheel = parsed.rr_vmotor / gearParams.gear_ratio;

fl_omega_wheel = fl_vwheel * (2*pi/60);  fr_omega_wheel = fr_vwheel * (2*pi/60);
rl_omega_wheel = rl_vwheel * (2*pi/60);  rr_omega_wheel = rr_vwheel * (2*pi/60);
fl_omega_motor = parsed.fl_vmotor * (2*pi/60);  fr_omega_motor = parsed.fr_vmotor * (2*pi/60);
rl_omega_motor = parsed.rl_vmotor * (2*pi/60);  rr_omega_motor = parsed.rr_vmotor * (2*pi/60);

% Clamp forces & mu-based brake torque, needed only to reconstruct Tbias_brake
front_piston_area = aeroParams.front_piston_count * pi * (aeroParams.front_piston_dia/2)^2;
rear_piston_area  = aeroParams.rear_piston_count  * pi * (aeroParams.rear_piston_dia/2)^2;
front_rotor_radius = aeroParams.front_rotor_dia / 2;
rear_rotor_radius  = aeroParams.rear_rotor_dia  / 2;

fl_clamp_force = (frontpressure * 6895) .* front_piston_area;
fr_clamp_force = (frontpressure * 6895) .* front_piston_area;
rl_clamp_force = (rearpressure  * 6895) .* rear_piston_area;
rr_clamp_force = (rearpressure  * 6895) .* rear_piston_area;

mu_front = interp1(aeroParams.mu_temp_table, aeroParams.mu_table, fr_temp_F, 'linear', 'extrap');
mu_front = max(min(mu_front, max(aeroParams.mu_table)), min(aeroParams.mu_table));
mu_rear  = interp1(aeroParams.mu_temp_table, aeroParams.mu_table, rr_temp_F, 'linear', 'extrap');
mu_rear  = max(min(mu_rear,  max(aeroParams.mu_table)), min(aeroParams.mu_table));

fl_Tbrake = -2 * mu_front .* fl_clamp_force .* front_rotor_radius;
fr_Tbrake = -2 * mu_front .* fr_clamp_force .* front_rotor_radius;
rl_Tbrake = -2 * mu_rear  .* rl_clamp_force .* rear_rotor_radius;
rr_Tbrake = -2 * mu_rear  .* rr_clamp_force .* rear_rotor_radius;
ftot_Tbrake = fl_Tbrake + fr_Tbrake;
rtot_Tbrake = rl_Tbrake + rr_Tbrake;
Tbias_brake = ftot_Tbrake ./ (rtot_Tbrake + ftot_Tbrake);
Tbias_brake(isnan(Tbias_brake)) = 0;

% Aero drag (DRS-dependent)
F_aero_open   = aeroParams.aero_open_a   * speed_mph.^2 + aeroParams.aero_open_b   * speed_mph + aeroParams.aero_open_c;
F_aero_closed = aeroParams.aero_closed_a * speed_mph.^2 + aeroParams.aero_closed_b * speed_mph + aeroParams.aero_closed_c;
F_aero = F_aero_open .* double(parsed.drs_state == 1) + F_aero_closed .* double(parsed.drs_state == 0);
F_aero = max(F_aero, 0);

% Regen power
decelerating_idx = parsed.accelx < 0;
fl_regen_power = min(fl_Tmotor .* fl_omega_motor, 0) .* double(decelerating_idx);
fr_regen_power = min(fr_Tmotor .* fr_omega_motor, 0) .* double(decelerating_idx);
rl_regen_power = min(rl_Tmotor .* rl_omega_motor, 0) .* double(decelerating_idx);
rr_regen_power = min(rr_Tmotor .* rr_omega_motor, 0) .* double(decelerating_idx);
total_regen_power = fl_regen_power + fr_regen_power + rl_regen_power + rr_regen_power;

derived.velx              = velx;
derived.frontpressure     = frontpressure;
derived.rearpressure      = rearpressure;
derived.fr_temp_F         = fr_temp_F;
derived.rr_temp_F         = rr_temp_F;
derived.Tbias_brake       = Tbias_brake;
derived.F_aero            = F_aero;
derived.total_regen_power = total_regen_power;
derived.fl_omega_wheel    = fl_omega_wheel;
derived.fr_omega_wheel    = fr_omega_wheel;
derived.rl_omega_wheel    = rl_omega_wheel;
derived.rr_omega_wheel    = rr_omega_wheel;
end


function Edrag = compute_edrag(t, velx, F_aero)
% Same as the "Calculate Aero Braking Energy" section of brake_temp_sim.m
d = zeros(size(t));
for n = 2:length(t)
    avgspd = (velx(n-1) + velx(n)) / 2;
    d(n)   = avgspd * (t(n) - t(n-1));
end
Edrag = zeros(size(t));
for k = 2:length(t)
    Edrag(k) = F_aero(k) * d(k);
end
end


function residual = brake_temp_residuals(p, padfrac_fun, datasets, ...
    VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
    I, WheelR, residualMode)
% Pools front + rear residuals across ALL datasets into one vector for
% lsqnonlin. x1,b1 (h_w) and the PadFrac params are shared/global across
% every dataset and both corners - only Tbias flips between front/rear.
%
% residualMode: 'all' (default) uses every sample. 'cooldown_only' zeros
% out residuals from steps that took the active-braking branch, so the
% objective only reflects passive-cooling behavior - used by the
% two-stage fit to identify h_w without PadFrac able to compensate for it.
if nargin < 11
    residualMode = 'all';
end

x1f_p = p(1); b1f_p = p(2); x1r_p = p(3); b1r_p = p(4); padfrac_params = p(5:end);

residual = [];
for k = 1:numel(datasets)
    ds = datasets(k);

    [predF_front, activeFront] = run_sim_opt(ds.t, ds.velx, ds.frontpressure, ds.fr_temp_F, ds.Tbias_brake, ...
        x1f_p, b1f_p, padfrac_fun, padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.fl_omega_wheel, ds.fr_omega_wheel, VehicleMass, RotorMass_front, RotorArea_front, I, WheelR, ds.TambK);

    [predF_rear, activeRear] = run_sim_opt(ds.t, ds.velx, ds.rearpressure, ds.rr_temp_F, 1 - ds.Tbias_brake, ...
        x1r_p, b1r_p, padfrac_fun, padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.rl_omega_wheel, ds.rr_omega_wheel, VehicleMass, RotorMass_rear, RotorArea_rear, I, WheelR, ds.TambK);

    residFront = predF_front(:) - ds.fr_temp_F(:);
    residRear  = predF_rear(:) - ds.rr_temp_F(:);

    if strcmp(residualMode, 'cooldown_only')
        residFront(activeFront) = 0;
        residRear(activeRear)   = 0;
    end

    residual = [residual; residFront; residRear]; %#ok<AGROW>
end

residual(~isfinite(residual)) = 0;  % safety net against a rare divergent step
end


function [RotorTempArrayF, activeMask] = run_sim_opt(t, velx, BrakePress, brakeTempArray, Tbias, ...
    x1_p, b1_p, padfrac_fun, padfrac_params, total_regen_power, Edrag, ...
    omega_wheel_L, omega_wheel_R, VehicleMass, RotorMass, RotorArea, I, WheelR, TambK) %#ok<INUSD>
% Same physics as run_sim() in brake_temp_sim.m, except PadFrac is now
% evaluated by an arbitrary function of (rotor temp [K], applied
% pressure [psi]) instead of a hardcoded linear-in-T formula.
%
% activeMask(i) is true where step i took the active-braking branch (i.e.
% PadFrac actually mattered for that step) and false where it took the
% passive-cooling branch. Used to isolate cooldown-only residuals for the
% two-stage fitting strategy below.

RotorTempArrayK = zeros(size(t));
RotorTempArrayK(1) = (brakeTempArray(1) - 32) * (5/9) + 273.15;
activeMask = false(size(t));
min_pressure = 5;  % psi, threshold to consider brakes applied

for i = 2:length(t)
    prevSpeed = velx(i-1);
    newSpeed  = velx(i);
    DS        = newSpeed - prevSpeed;
    prevTemp  = RotorTempArrayK(i-1);
    tbrake    = t(i) - t(i-1);
    if tbrake <= 0
        RotorTempArrayK(i) = prevTemp;
        continue
    end

    h_w      = x1_p * velx(i) + b1_p;
    Rrotor   = 1 / (h_w * RotorArea);
    PadFrac  = max(min(padfrac_fun(prevTemp, BrakePress(i), padfrac_params), 1), 0);
    SpecHeat = (0.0005 * prevTemp + 0.2813) * 1000;

    if DS < -2.5
        DS = 0;  % skip - exceeds 2.5g deceleration limit
    end

    if DS < 0 && BrakePress(i) > min_pressure
        activeMask(i) = true;
        Energy1 = 0.5 * VehicleMass * (prevSpeed^2 - newSpeed^2);
        omegaP  = (omega_wheel_L(i-1) + omega_wheel_R(i-1)) / 2;
        omegaN  = (omega_wheel_L(i)   + omega_wheel_R(i))   / 2;
        Energy2 = 4 * (0.5 * I * (omegaP^2 - omegaN^2));
        Energy  = Energy1 + Energy2;

        regen_energy    = abs(total_regen_power(i)) * tbrake;
        friction_energy = max(Energy - regen_energy, 0);
        AeroFrac  = min(Edrag(i) / max(Energy, 1), 1);
        BrakeFrac = 0.82;

        CorrectedEnergy    = friction_energy * 0.5 * Tbias(i) * (1 - AeroFrac) * (1 - PadFrac) * BrakeFrac * 0.8;
        deltaTK             = CorrectedEnergy / (RotorMass * SpecHeat);
        RotorTempArrayK(i)  = deltaTK + prevTemp;

        qout       = (RotorTempArrayK(i) - TambK) / Rrotor;
        Eout       = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = RotorTempArrayK(i) - deltaTKout;
    else
        qout       = (prevTemp - TambK) / Rrotor;
        Eout       = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = prevTemp - deltaTKout;
    end
end

RotorTempArrayF = ((RotorTempArrayK - 273.15) * (9/5)) + 32;
end