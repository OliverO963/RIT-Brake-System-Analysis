%% ================================================================
%  Rotor Temperature Check  (single rotor mass / area, no optimization)
% ================================================================
%  Evaluates ONE rotor mass and surface area per axle - RotorMass_front/rear
%  and RotorArea_front/rear from the Parameters sheet of VehicleParameters.xlsx
%  - so the temperature performance of a specific design can be checked.
%
%  For every previous-vehicle dataset (same import as BrakeCoeffOptimizer.m)
%  one figure per dataset, front and rear panels, three curves:
%    1. Recorded rotor temperature (old rotors, old vehicle)
%    2. Predicted: NEW rotor, PREVIOUS vehicle's regen   (run_sim_opt, unchanged)
%    3. Predicted: NEW rotor, NEW vehicle parameters and regen characteristics
%       derived from the vehicle-dynamics simulation CSV, applied to that
%       dataset's drive cycle.
%  One extra figure: NEW rotor on the simulated endurance cycle (the CSV),
%  driven directly by its FrictionBrakeForce_N column.
%
%  Established coefficients (h_w, PadFrac) are READ from the sheet, never fit.
%  All masses/areas/energies are PER SINGLE ROTOR (one corner).
% ================================================================

clc; close all

%% ================== CONFIG ==================
cfg.paramFile = 'VehicleParameters.xlsx';
cfg.simCsv    = 'F35_8HardLapsSim.csv';   % vehicle-dynamics endurance simulation
cfg.files     = {};                       % {} -> file picker; or full paths to .txt datasets
cfg.csvSpeedUnits = 'mph';                % Speed_mps column is labelled m/s but is MPH (see check printed below)
cfg.makePlots = true;
cfg.showDiagnostics = true;               % extra figure: derived regen curve and aero fit vs the CSV

% Optional local overrides (e.g. paths) - runs if rotor_temp_overrides.m is on the path.
if exist('rotor_temp_overrides', 'file') == 2
    rotor_temp_overrides;
end

%% ================== PARAMETERS (from the sheet) ==================
if exist(cfg.paramFile, 'file') ~= 2
    [f, p] = uigetfile('*.xlsx', 'Select VehicleParameters.xlsx');
    if isequal(f, 0), error('No parameter file selected.'); end
    cfg.paramFile = fullfile(p, f);
end
P = read_params(cfg.paramFile);

mF = P.RotorMass_front;  AF = P.RotorArea_front;     % NEW rotor under test, ONE front rotor
mR = P.RotorMass_rear;   AR = P.RotorArea_rear;      % NEW rotor under test, ONE rear rotor
x1f = P.hw_x1_front;  b1f = P.hw_b1_front;           % h_w = x1*v + b1  [W/m^2-K], v in m/s
x1r = P.hw_x1_rear;   b1r = P.hw_b1_rear;
padfrac_fun    = @(T, Pr, p) p(1) .* T + p(2);       % PadFrac(T_K); sheet model is pressure-independent
padfrac_params = [P.padfrac_x2, P.padfrac_b2];

fprintf('Rotor under test: front %.4f kg / %.5f m^2 | rear %.4f kg / %.5f m^2\n', mF, AF, mR, AR);
fprintf('h_w front %.4f*v + %.4f | rear %.4f*v + %.4f | PadFrac = %.6f*T_K + %.4f\n', ...
    x1f, b1f, x1r, b1r, P.padfrac_x2, P.padfrac_b2);
fprintf('Target %.0f F | new-vehicle Tbias_front %.2f | m %.1f kg | I_corner %.3f kg*m^2 | WheelR %.3f m\n', ...
    P.T_target_F, P.Tbias_front, P.VehicleMass, P.I_corner, P.WheelR);

%% ================== NEW-VEHICLE SIMULATION CSV ==================
if exist(cfg.simCsv, 'file') ~= 2
    [f, p] = uigetfile('*.csv', 'Select the vehicle-dynamics simulation CSV');
    if isequal(f, 0), error('No simulation CSV selected.'); end
    cfg.simCsv = fullfile(p, f);
end
fid = fopen(cfg.simCsv);  hdr = strsplit(strtrim(fgetl(fid)), ',');  fclose(fid);
M = readmatrix(cfg.simCsv, 'NumHeaderLines', 1);
col = @(name) get_col(M, hdr, name);
S.t      = col('Time_s');
S.lap    = col('Lap');
S.rpm    = col('MotorRPM');
S.avail  = col('AvailableRegenMagnitude_N');
S.req    = col('RequiredBrakingMagnitude_N');
S.fric   = col('FrictionBrakeForce_N');
S.aero   = col('AeroDragMagnitude_N');
spd = col('Speed_mps');
if strcmp(cfg.csvSpeedUnits, 'mph'), S.v = spd * 0.44704; else, S.v = spd; end

% Unit check: wheel radius implied by MotorRPM and speed should be near WheelR.
om = S.rpm / P.gear_ratio * 2*pi/60;  ok = om > 20;
Rimp = median(S.v(ok) ./ om(ok));
fprintf('\nCSV: %d rows, %.1f s, %d laps | implied wheel radius %.4f m vs sheet %.3f m (%+.1f%%)\n', ...
    numel(S.t), S.t(end) - S.t(1), max(S.lap), Rimp, P.WheelR, 100 * (Rimp / P.WheelR - 1));
if abs(Rimp / P.WheelR - 1) > 0.15
    warning('Implied wheel radius is >15%% off WheelR: check cfg.csvSpeedUnits (CSV speed units).');
end

% Characteristics derived from the CSV:
%   aero drag      F = aero_k * v^2   (least squares, v in m/s)
%   regen available F_avail(motor RPM)  (single-valued and monotone in the CSV)
aero_k = sum(S.aero .* S.v.^2) / sum(S.v.^4);
aero_rms = sqrt(mean((S.aero - aero_k * S.v.^2).^2));
[rpmU, availU] = regen_curve(S.rpm, S.avail);
fprintf('Derived aero: F = %.4f*v^2 N (fit rms %.2f N, i.e. rho*CdA = %.3f kg/m)\n', aero_k, aero_rms, 2 * aero_k);
fprintf('Derived regen: available %.0f N at %.0f RPM down to %.0f N at %.0f RPM (%d lookup points)\n', ...
    availU(1), rpmU(1), availU(end), rpmU(end), numel(rpmU));

%% ================== LOAD PREVIOUS-VEHICLE DATASETS ==================
if isempty(cfg.files)
    [files, path] = uigetfile('*.txt', 'Select dataset file(s)', 'MultiSelect', 'on');
    if isequal(files, 0), error('No files selected.'); end
    if ischar(files), files = {files}; end
    fullfiles = fullfile(path, files);
else
    fullfiles = cfg.files;
end

datasets = struct('name', {}, 't', {}, 'velx', {}, 'frontpressure', {}, 'rearpressure', {}, ...
    'fr_temp_F', {}, 'rr_temp_F', {}, 'Tbias_brake', {}, 'Edrag', {}, ...
    'front_axle_regen_power', {}, 'rear_axle_regen_power', {}, 'omega_all', {}, 'TambK', {});

for k = 1:numel(fullfiles)
    [~, shortName, ext] = fileparts(fullfiles{k});
    fprintf('\nLoading dataset %d/%d: %s%s\n', k, numel(fullfiles), shortName, ext);
    raw = readmatrix(fullfiles{k}, 'FileType', 'text', 'Delimiter', '\t', ...
        'NumHeaderLines', 6, 'TreatAsMissing', {'', 'NaN'});
    switch size(raw, 2)
        case 27, fmt = 'drs_imu';
        case 22, fmt = 'drs';
        case 21, fmt = 'no_drs';
        otherwise
            error('Unrecognized column count (%d) in %s. Expected 21, 22, or 27.', size(raw, 2), fullfiles{k});
    end
    if size(raw, 1) < 3, warning('  %s has < 3 samples - skipping.', shortName); continue; end

    parsed  = parse_dataset_columns(raw, fmt);
    derived = compute_derived_quantities(parsed, P);
    Edrag   = compute_edrag(parsed.t, derived.velx, derived.F_aero);

    firstF = [derived.fr_temp_F(1), derived.rr_temp_F(1)];
    firstF = firstF(isfinite(firstF));
    if isempty(firstF), TambC_ds = P.TambC_fallback; else, TambC_ds = (mean(firstF) - 32) * 5/9; end

    idx = numel(datasets) + 1;
    datasets(idx).name = [shortName ext];
    datasets(idx).t = parsed.t;                     datasets(idx).velx = derived.velx;
    datasets(idx).frontpressure = derived.frontpressure;  datasets(idx).rearpressure = derived.rearpressure;
    datasets(idx).fr_temp_F = derived.fr_temp_F;    datasets(idx).rr_temp_F = derived.rr_temp_F;
    datasets(idx).Tbias_brake = derived.Tbias_brake;      datasets(idx).Edrag = Edrag;
    datasets(idx).front_axle_regen_power = derived.front_axle_regen_power;
    datasets(idx).rear_axle_regen_power  = derived.rear_axle_regen_power;
    datasets(idx).omega_all = [derived.fl_omega_wheel(:), derived.fr_omega_wheel(:), ...
                               derived.rl_omega_wheel(:), derived.rr_omega_wheel(:)];
    datasets(idx).TambK = TambC_ds + 273.15;
    fprintf('  %s | %d samples, %.1f s | ambient %.1f F\n', fmt, numel(parsed.t), ...
        parsed.t(end) - parsed.t(1), TambC_ds * 9/5 + 32);
end
nDS = numel(datasets);
if nDS == 0, error('No usable datasets loaded.'); end

%% ================== SIMULATE: previous-vehicle datasets ==================
res = struct('b_front', {}, 'b_rear', {}, 'c_front', {}, 'c_rear', {}, 'regen_share', {});
for k = 1:nDS
    ds = datasets(k);
    % (2) NEW rotor, PREVIOUS vehicle's regen: run_sim_opt exactly as in BrakeCoeffOptimizer.m
    res(k).b_front = run_sim_opt(ds.t, ds.velx, ds.frontpressure, ds.fr_temp_F, ds.Tbias_brake, ...
        x1f, b1f, padfrac_fun, padfrac_params, ds.front_axle_regen_power, ds.Edrag, ds.omega_all, ...
        P.VehicleMass, mF, AF, P.I_corner, ds.TambK, P);
    res(k).b_rear = run_sim_opt(ds.t, ds.velx, ds.rearpressure, ds.rr_temp_F, 1 - ds.Tbias_brake, ...
        x1r, b1r, padfrac_fun, padfrac_params, ds.rear_axle_regen_power, ds.Edrag, ds.omega_all, ...
        P.VehicleMass, mR, AR, P.I_corner, ds.TambK, P);
    % (3) NEW rotor, NEW vehicle + CSV-derived regen, same drive cycle
    [Ef, Er, res(k).regen_share] = new_vehicle_friction(ds, P, aero_k, rpmU, availU);
    [res(k).c_front, Rrotor_cf] = run_sim_force(ds.t, ds.velx, Ef, ds.fr_temp_F(1), ds.TambK, x1f, b1f, ...
        padfrac_fun, padfrac_params, mF, AF, P);
    [res(k).c_rear, Rrotor_cr]  = run_sim_force(ds.t, ds.velx, Er, ds.rr_temp_F(1), ds.TambK, x1r, b1r, ...
        padfrac_fun, padfrac_params, mR, AR, P);
end

%% ================== SIMULATE: endurance-pace CSV ==================
[Ef_csv, Er_csv, csv_regen_share] = csv_friction(S, P);
T0 = P.T_init_F;  Tamb_csv = P.TambC_fallback + 273.15;
[Tcsv_front, Rrotor_csvf] = run_sim_force(S.t, S.v, Ef_csv, T0, Tamb_csv, x1f, b1f, padfrac_fun, padfrac_params, mF, AF, P);
[Tcsv_rear, Rrotor_csvr]  = run_sim_force(S.t, S.v, Er_csv, T0, Tamb_csv, x1r, b1r, padfrac_fun, padfrac_params, mR, AR, P);

%% ================== REPORT ==================
fprintf('\n================ PEAK ROTOR TEMPERATURE (deg F) ================\n');
fprintf('%-44s %-6s %10s %14s %14s %12s\n', 'Dataset', 'Axle', 'Recorded', 'New+old regen', 'New+new regen', 'New regen %');
for k = 1:nDS
    fprintf('%-44s %-6s %10.1f %14.1f %14.1f %12.0f\n', datasets(k).name, 'front', max(datasets(k).fr_temp_F), ...
        max(res(k).b_front), max(res(k).c_front), 100 * res(k).regen_share);
    fprintf('%-44s %-6s %10.1f %14.1f %14.1f %12s\n', '', 'rear', max(datasets(k).rr_temp_F), ...
        max(res(k).b_rear), max(res(k).c_rear), '');
end
fprintf('\n================ SIMULATED ENDURANCE (%s) ================\n', cfg.simCsv);
fprintf('Target %.0f F. CSV regen share of required braking energy: %.0f%%\n', P.T_target_F, 100 * csv_regen_share);
fprintf('front: peak %.1f F (%+.1f vs target), %.1f s above target\n', max(Tcsv_front), max(Tcsv_front) - P.T_target_F, ...
    sum(diff(S.t(:)) .* (Tcsv_front(2:end) > P.T_target_F)));
fprintf('rear : peak %.1f F (%+.1f vs target), %.1f s above target\n', max(Tcsv_rear), max(Tcsv_rear) - P.T_target_F, ...
    sum(diff(S.t(:)) .* (Tcsv_rear(2:end) > P.T_target_F)));

%% ================== PLOTS ==================
if cfg.makePlots
    for k = 1:nDS
        ds = datasets(k);
        figure('Name', ds.name);
        for a = 1:2
            if a == 1
                meas = ds.fr_temp_F; b = res(k).b_front; c = res(k).c_front; nm = 'Front';
            else
                meas = ds.rr_temp_F; b = res(k).b_rear;  c = res(k).c_rear;  nm = 'Rear';
            end
            subplot(2, 1, a);  hold on;  grid on;
            plot(ds.t, meas, 'k-', 'DisplayName', 'Recorded (old rotor)');
            plot(ds.t, b, 'b-', 'DisplayName', 'New rotor, previous-vehicle regen');
            plot(ds.t, c, 'r-', 'DisplayName', 'New rotor, new-vehicle parameters + regen');
            yline(P.T_target_F, ':', 'Target', 'Color', [0.4 0.4 0.4]);
            xlabel('Time (s)');  ylabel([nm ' rotor temp (deg F)']);
            legend('Location', 'best');  title([nm ' rotor']);
        end
        sgtitle(ds.name, 'Interpreter', 'none');
    end

    figure('Name', 'Simulated endurance - new rotor');
    subplot(2, 1, 1);  hold on;  grid on;
    plot(S.t, Tcsv_front, 'r-', 'DisplayName', sprintf('Front (%.3f kg / %.4f m^2)', mF, AF));
    % yline(P.T_target_F, ':', 'Target', 'Color', [0.4 0.4 0.4]);
    xlabel('Time (s)');  ylabel('Front rotor temp (deg F)');  legend('Location', 'best');
    title(sprintf('Front: peak %.0f F', max(Tcsv_front)));
    subplot(2, 1, 2);  hold on;  grid on;
    plot(S.t, Tcsv_rear, 'b-', 'DisplayName', sprintf('Rear (%.3f kg / %.4f m^2)', mR, AR));
    % yline(P.T_target_F, ':', 'Target', 'Color', [0.4 0.4 0.4]);
    xlabel('Time (s)');  ylabel('Rear rotor temp (deg F)');  legend('Location', 'best');
    title(sprintf('Rear: peak %.0f F', max(Tcsv_rear)));
    sgtitle('Simulated endurance pace (CSV), new rotors');

    if cfg.showDiagnostics
        figure('Name', 'Derived characteristics vs CSV');
        subplot(1, 2, 1);  hold on;  grid on;
        plot(S.rpm, S.avail, 'k.', 'DisplayName', 'CSV AvailableRegenMagnitude');
        plot(rpmU, availU, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Derived F_{avail}(RPM)');
        xlabel('Motor RPM');  ylabel('Available regen force (N)');  legend('Location', 'best');
        subplot(1, 2, 2);  hold on;  grid on;
        plot(S.v, S.aero, 'k.', 'DisplayName', 'CSV AeroDragMagnitude');
        vv = linspace(0, max(S.v), 100);
        plot(vv, aero_k * vv.^2, 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('%.3f v^2', aero_k));
        xlabel('Speed (m/s)');  ylabel('Aero drag (N)');  legend('Location', 'best');
    end
end


%% ================================================================
%  LOCAL FUNCTIONS
% ================================================================

function P = read_params(file)
% Reads Name/Value pairs from the Parameters sheet and builds the same
% fields the data-import functions expect. Piston area per side: an
% opposed caliper's far side reacts the load, so it uses half the count.
c = readcell(file, 'Sheet', 'Parameters');
P = struct();
for r = 2:size(c, 1)
    nm = c{r, 1};  val = c{r, 2};
    if (ischar(nm) || isstring(nm)) && isnumeric(val) && isscalar(val)
        P.(char(nm)) = val;
    end
end
mu = readmatrix(file, 'Sheet', 'MuTable', 'NumHeaderLines', 1);
P.mu_temp_table = mu(:, 1)';  P.mu_table = mu(:, 2)';
P.front_piston_area_per_side = P.front_piston_count / (1 + P.front_caliper_opposed) * pi * (P.front_piston_dia/2)^2;
P.rear_piston_area_per_side  = P.rear_piston_count  / (1 + P.rear_caliper_opposed)  * pi * (P.rear_piston_dia/2)^2;
req = {'VehicleMass','I_corner','WheelR','gear_ratio','TambC_fallback','RotorMass_front','RotorMass_rear', ...
    'RotorArea_front','RotorArea_rear','BrakeFrac','CalibrationFactor','min_pressure','T_target_F','T_init_F', ...
    'Tbias_front','hw_x1_front','hw_b1_front','hw_x1_rear','hw_b1_rear','padfrac_x2','padfrac_b2', ...
    'press_front_slope','press_front_offset','press_rear_slope','press_rear_offset','temp_adc_slope', ...
    'temp_adc_offset','velx_threshold','front_pad_mean_radius','rear_pad_mean_radius', ...
    'aero_open_a','aero_open_b','aero_open_c','aero_closed_a','aero_closed_b','aero_closed_c'};
missing = req(~isfield(P, req));
if ~isempty(missing)
    error('Parameters sheet is missing: %s', strjoin(missing, ', '));
end
end

function v = get_col(M, hdr, name)
j = find(strcmp(hdr, name), 1);
if isempty(j), error('CSV column "%s" not found.', name); end
v = M(:, j);
end

function [rpmU, availU] = regen_curve(rpm, avail)
% Available regen force as a function of motor RPM, from the CSV (the CSV's
% availability is single-valued and monotone-decreasing in RPM, not in vehicle speed).
key = round(rpm * 10) / 10;
[rpmU, ~, j] = unique(key);
availU = accumarray(j, avail, [], @mean);
rpmU = rpmU(:);  availU = availU(:);
end

function F = regen_avail(rpm, rpmU, availU)
% Lookup, clamped to the CSV's RPM range (no extrapolation beyond it).
F = interp1(rpmU, availU, min(max(rpm, rpmU(1)), rpmU(end)), 'linear');
end

function [Ef, Er, regen_share] = new_vehicle_friction(ds, P, aero_k, rpmU, availU)
% NEW vehicle on a previous-vehicle drive cycle. Per decelerating step:
%   required braking force  F_req = m_eq*decel - aero      (m_eq = m + 4*I_corner/R^2)
%   regen force             F_regen = min(F_req, F_avail(motor RPM))  [CSV: applied = min(required, available)]
%   friction force          F_fric = F_req - F_regen
%   vehicle friction energy = F_fric * distance; corner energy = x Tbias_front (rear: 1-Tbias) x 0.5
% A step counts only if the previous vehicle's driver was braking (friction pressure or regen present),
% which keeps speed-sensor noise from creating phantom braking.
t = ds.t(:);  v = ds.velx(:);
dt = [0; diff(t)];  dv = [0; diff(v)];
vavg = [0; (v(1:end-1) + v(2:end)) / 2];
decel = -dv ./ max(dt, eps);
valid = dt > 0 & dv < 0 & dv >= -2.5;           % same glitch rule as run_sim_opt (DS < -2.5 skipped)
demand = (ds.frontpressure(:) > P.min_pressure) | (ds.rearpressure(:) > P.min_pressure) | ...
         (ds.front_axle_regen_power(:) + ds.rear_axle_regen_power(:)) < 0;
m_eq = P.VehicleMass + 4 * P.I_corner / P.WheelR^2;
F_req = max(m_eq * decel - aero_k * vavg.^2, 0) .* (valid & demand);
rpm = mean(abs(ds.omega_all), 2) * P.gear_ratio * 60 / (2*pi);
F_regen = min(F_req, regen_avail(rpm, rpmU, availU));
F_fric = F_req - F_regen;
dist = vavg .* dt;
Ev = F_fric .* dist;
Ef = Ev * P.Tbias_front * 0.5;
Er = Ev * (1 - P.Tbias_front) * 0.5;
regen_share = sum(F_regen .* dist) / max(sum(F_req .* dist), eps);
end

function [Ef, Er, regen_share] = csv_friction(S, P)
% Endurance-pace CSV: friction force is given directly (FrictionBrakeForce_N, a contact-patch force).
dt = [0; diff(S.t(:))];
vavg = [0; (S.v(1:end-1) + S.v(2:end)) / 2];
Ev = abs(S.fric(:)) .* vavg .* dt;
Ef = Ev * P.Tbias_front * 0.5;
Er = Ev * (1 - P.Tbias_front) * 0.5;
regen_share = 1 - sum(abs(S.fric(:)) .* vavg .* dt) / max(sum(S.req(:) .* vavg .* dt), eps);
end

function [T_F, Rrotor] = run_sim_force(t, v, E_corner, T0_F, TambK, x1, b1, padfrac_fun, padfrac_params, m, A, P)
% Rotor thermal update IDENTICAL to run_sim_opt (heat in, then convection out
% using the updated temperature), except the heat input is supplied directly as
% per-step friction energy at ONE corner (J) instead of being rebuilt from
% KE loss. Passive cooling when E_corner(i) = 0.
n = numel(t);
TK = zeros(n, 1);  TK(1) = (T0_F - 32) * 5/9 + 273.15;
Rrotor = zeros(n, 1); Rrotor(1) = 1/ ((x1 * v(1) + b1) * A);
for i = 2:n
    tb = t(i) - t(i-1);  prev = TK(i-1);
    if tb <= 0, TK(i) = prev; continue; end
    Rrotor(i)  = 1 / ((x1 * v(i) + b1) * A);
    PadFrac = max(min(padfrac_fun(prev, 0, padfrac_params), 1), 0);
    cp = (0.0005 * prev + 0.2813) * 1000;
    if E_corner(i) > 0
        TK(i) = prev + E_corner(i) * (1 - PadFrac) * P.BrakeFrac * P.CalibrationFactor / (m * cp);
        TK(i) = TK(i) - ((TK(i) - TambK) / Rrotor(i)) * tb / (m * cp);
    else
        TK(i) = prev - ((prev - TambK) / Rrotor(i)) * tb / (m * cp);
    end
end
T_F = (TK - 273.15) * 9/5 + 32;
end

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

function derived = compute_derived_quantities(parsed, VP)
% Reconstructs everything brake_temp_sim needs: cleaned velocity,
% front/rear pressure (psi), front/rear rotor temp (degF), Tbias_brake,
% aero drag force, per-axle regen power, and per-corner wheel omega.
% P is the centralized parameter struct from loadVehicleParams().

velx = parsed.velx;
velx(abs(velx) > VP.velx_threshold) = 0;
velx(isnan(velx)) = 0;
velx(velx < 0)    = 0;
speed_mph = velx * 2.23694;

frontpressure = max(VP.press_front_slope * parsed.frontpressure_adc - VP.press_front_offset, 0);  % psi
rearpressure  = max(VP.press_rear_slope  * parsed.rearpressure_adc  - VP.press_rear_offset,  0);  % psi

fr_temp_C = VP.temp_adc_slope * (parsed.fr_temp_adc - VP.temp_adc_offset);
rr_temp_C = VP.temp_adc_slope * (parsed.rr_temp_adc - VP.temp_adc_offset);
fr_temp_F = fr_temp_C * (9/5) + 32;
rr_temp_F = rr_temp_C * (9/5) + 32;

fl_Tmotor = parsed.fl_Tmotor_Mn / 100 * 9.8;  fr_Tmotor = parsed.fr_Tmotor_Mn / 100 * 9.8;
rl_Tmotor = parsed.rl_Tmotor_Mn / 100 * 9.8;  rr_Tmotor = parsed.rr_Tmotor_Mn / 100 * 9.8;

fl_vwheel = parsed.fl_vmotor / VP.gear_ratio;  fr_vwheel = parsed.fr_vmotor / VP.gear_ratio;
rl_vwheel = parsed.rl_vmotor / VP.gear_ratio;  rr_vwheel = parsed.rr_vmotor / VP.gear_ratio;

fl_omega_wheel = fl_vwheel * (2*pi/60);  fr_omega_wheel = fr_vwheel * (2*pi/60);
rl_omega_wheel = rl_vwheel * (2*pi/60);  rr_omega_wheel = rr_vwheel * (2*pi/60);
fl_omega_motor = parsed.fl_vmotor * (2*pi/60);  fr_omega_motor = parsed.fr_vmotor * (2*pi/60);
rl_omega_motor = parsed.rl_vmotor * (2*pi/60);  rr_omega_motor = parsed.rr_vmotor * (2*pi/60);

% Clamp forces & mu-based brake torque, needed only to reconstruct
% Tbias_brake. Clamp force uses the SINGLE-SIDE piston area (loadVehicleParams
% halves the per-caliper piston count for opposed calipers) - the caliper is
% a closed force loop, so the far side reacts the load rather than adding to
% it. The factor of 2 below is the two friction FACES, not the two sides'
% pistons. Lever arm is the documented MEAN PAD RADIUS, not the rotor OD.
% NOTE: both of these cancel out of Tbias_brake (it is a ratio and every
% term carries the same factor), so they do not change this script's fit -
% they matter for absolute torque, which padSpecificPower.m consumes.
front_piston_area_side = VP.front_piston_area_per_side;
rear_piston_area_side  = VP.rear_piston_area_per_side;
front_pad_radius = VP.front_pad_mean_radius;
rear_pad_radius  = VP.rear_pad_mean_radius;

fl_clamp_force = (frontpressure * 6895) .* front_piston_area_side;
fr_clamp_force = (frontpressure * 6895) .* front_piston_area_side;
rl_clamp_force = (rearpressure  * 6895) .* rear_piston_area_side;
rr_clamp_force = (rearpressure  * 6895) .* rear_piston_area_side;

mu_front = interp1(VP.mu_temp_table, VP.mu_table, fr_temp_F, 'linear', 'extrap');
mu_front = max(min(mu_front, max(VP.mu_table)), min(VP.mu_table));
mu_rear  = interp1(VP.mu_temp_table, VP.mu_table, rr_temp_F, 'linear', 'extrap');
mu_rear  = max(min(mu_rear,  max(VP.mu_table)), min(VP.mu_table));

fl_Tbrake = -2 * mu_front .* fl_clamp_force .* front_pad_radius;
fr_Tbrake = -2 * mu_front .* fr_clamp_force .* front_pad_radius;
rl_Tbrake = -2 * mu_rear  .* rl_clamp_force .* rear_pad_radius;
rr_Tbrake = -2 * mu_rear  .* rr_clamp_force .* rear_pad_radius;
ftot_Tbrake = fl_Tbrake + fr_Tbrake;
rtot_Tbrake = rl_Tbrake + rr_Tbrake;
Tbias_brake = ftot_Tbrake ./ (rtot_Tbrake + ftot_Tbrake);
Tbias_brake(isnan(Tbias_brake)) = 0;

% Aero drag (DRS-dependent)
F_aero_open   = VP.aero_open_a   * speed_mph.^2 + VP.aero_open_b   * speed_mph + VP.aero_open_c;
F_aero_closed = VP.aero_closed_a * speed_mph.^2 + VP.aero_closed_b * speed_mph + VP.aero_closed_c;
F_aero = F_aero_open .* double(parsed.drs_state == 1) + F_aero_closed .* double(parsed.drs_state == 0);
F_aero = max(F_aero, 0);

% Regen power
decelerating_idx = parsed.accelx < 0;
fl_regen_power = min(fl_Tmotor .* fl_omega_motor, 0) .* double(decelerating_idx);
fr_regen_power = min(fr_Tmotor .* fr_omega_motor, 0) .* double(decelerating_idx);
rl_regen_power = min(rl_Tmotor .* rl_omega_motor, 0) .* double(decelerating_idx);
rr_regen_power = min(rr_Tmotor .* rr_omega_motor, 0) .* double(decelerating_idx);

derived.velx              = velx;
derived.frontpressure     = frontpressure;
derived.rearpressure      = rearpressure;
derived.fr_temp_F         = fr_temp_F;
derived.rr_temp_F         = rr_temp_F;
derived.Tbias_brake       = Tbias_brake;
derived.F_aero            = F_aero;
% Regen is measured PER CORNER, so it is kept per axle rather than pooled to
% a single vehicle-level total. Pooling it and subtracting before the
% Tbias/corner split let front regen offset REAR rotor heat in proportion to
% brake bias, which is not physical.
derived.front_axle_regen_power = fl_regen_power + fr_regen_power;
derived.rear_axle_regen_power  = rl_regen_power + rr_regen_power;
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

function [RotorTempArrayF, activeMask] = run_sim_opt(t, velx, BrakePress, brakeTempArray, Tbias, ...
    x1_p, b1_p, padfrac_fun, padfrac_params, axle_regen_power, Edrag, ...
    omega_all, VehicleMass, RotorMass, RotorArea, I, TambK, VP)
% Same physics as run_sim() in brake_temp_sim.m, except PadFrac is now
% evaluated by an arbitrary function of (rotor temp [K], applied
% pressure [psi]) instead of a hardcoded linear-in-T formula.
%
% BASIS: every sizing input is per single component. The energy chain is
%   whole vehicle -> x Tbias -> one axle -> x0.5 -> ONE CORNER = ONE ROTOR,
% so the corner energy divides by one RotorMass directly.
%
% `omega_all` is a 4-column [fl fr rl rr] matrix of per-corner wheel speeds,
% used to sum rotational KE strictly per corner. `axle_regen_power` is THIS
% axle's regen only, so regen is attributed to the corners that produced it.
%
% activeMask(i) is true where step i took the active-braking branch (i.e.
% PadFrac actually mattered for that step) and false where it took the
% passive-cooling branch. Used to isolate cooldown-only residuals for the
% two-stage fitting strategy below.

RotorTempArrayK = zeros(size(t));
RotorTempArrayK(1) = (brakeTempArray(1) - 32) * (5/9) + 273.15;
activeMask = false(size(t));
min_pressure = VP.min_pressure;   % psi, threshold to consider brakes applied
corner_split = 0.5;              % one axle -> ONE corner (= one rotor)

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
        % Rotational KE summed strictly PER CORNER: I is one corner's
        % inertia, so each corner contributes using its OWN omega.
        Energy2 = sum(0.5 * I * (omega_all(i-1,:).^2 - omega_all(i,:).^2));
        Energy  = Energy1 + Energy2;

        AeroFrac  = min(Edrag(i) / max(Energy, 1), 1);

        % Split the vehicle's KE loss down to THIS CORNER first, then remove
        % aero and this corner's own regen (previously regen was pooled
        % across all four corners and subtracted before the split).
        Energy_corner   = Energy * Tbias(i) * corner_split * (1 - AeroFrac);
        regen_energy    = abs(axle_regen_power(i)) * tbrake * corner_split;
        friction_energy = max(Energy_corner - regen_energy, 0);

        CorrectedEnergy    = friction_energy * (1 - PadFrac) * VP.BrakeFrac * VP.CalibrationFactor;
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

