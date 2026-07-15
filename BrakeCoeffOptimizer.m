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
TambC           = 25;       % degC

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
    'rl_omega_wheel', {}, 'rr_omega_wheel', {});

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
    prompt = {sprintf('Start time (s)  [data range %.1f - %.1f]:', t_full(1), t_full(end)), 'End time (s):'};
    defaultAns = {num2str(t_full(1)), num2str(t_full(end))};
    answer = inputdlg(prompt, sprintf('Time crop: %s', files{k}), 1, defaultAns);
    if isempty(answer)
        t_start_k = t_full(1); t_end_k = t_full(end);
    else
        t_start_k = str2double(answer{1});
        t_end_k   = str2double(answer{2});
    end
    raw = raw(t_full >= t_start_k & t_full <= t_end_k, :);
    if size(raw, 1) < 3
        warning('  Dataset %s has < 3 samples after crop - skipping.', files{k});
        continue
    end

    parsed  = parse_dataset_columns(raw, fmt);
    derived = compute_derived_quantities(parsed, gearParams, aeroParams);
    Edrag   = compute_edrag(parsed.t, derived.velx, derived.F_aero);

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
models(2).lb   = [-x2_bound, -x3_bound, -b2_bound];
models(2).ub   = [ x2_bound,  x3_bound,  b2_bound];
models(2).x0   = [x2_seed, 0, b2_seed];

models(3).name = 'Linear with T*P interaction';
models(3).fun  = @(T, P, p) p(1).*T + p(2).*P + p(3).*T.*P + p(4);
models(3).lb   = [-x2_bound, -x3_bound, -x4_bound, -b2_bound];
models(3).ub   = [ x2_bound,  x3_bound,  x4_bound,  b2_bound];
models(3).x0   = [x2_seed, 0, 0, b2_seed];

models(4).name = 'Quadratic in T, linear in P';
quad_bound = x2_bound / dT;
models(4).fun  = @(T, P, p) p(1).*T + p(2).*T.^2 + p(3).*P + p(4);
models(4).lb   = [-x2_bound, -quad_bound, -x3_bound, -b2_bound];
models(4).ub   = [ x2_bound,  quad_bound,  x3_bound,  b2_bound];
models(4).x0   = [x2_seed, 0, 0, b2_seed];

%% ================== RUN OPTIMIZATION FOR EACH MODEL ==================
opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxFunctionEvaluations', 1000, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

results = struct('name', {}, 'params', {}, 'h_w', {}, 'padfrac_params', {}, ...
    'rmse_F', {}, 'avg_pct_err', {}, 'sse', {}, 'nresid', {}, 'nparams', {}, 'aicc', {});

TambK = TambC + 273.15;

for m = 1:numel(models)
    fprintf('\n=== Fitting model %d/%d: %s ===\n', m, numel(models), models(m).name);

    x0 = [x1_seed, b1_seed, x1_seed, b1_seed, models(m).x0];
    lb = [x1_lb, b1_lb, x1_lb, b1_lb, models(m).lb];
    ub = [x1_ub, b1_ub, x1_lb, b1_lb, models(m).ub];

    resFun = @(p) brake_temp_residuals(p, models(m).fun, datasets, ...
        VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
        I, WheelR, TambK);

    [p_opt, resnorm, residual] = lsqnonlin(resFun, x0, lb, ub, opts);

    residual = residual(isfinite(residual));
    nResid   = numel(residual);
    sse      = sum(residual.^2);
    rmse     = sqrt(sse / nResid);
    kParams  = numel(p_opt);
    % AICc: lower is better; penalizes extra parameters so a model isn't
    % favored just because it has more knobs to turn.
    aicc = nResid*log(sse/nResid) + 2*kParams + (2*kParams*(kParams+1)) / max(nResid - kParams - 1, 1);

    results(m).name           = models(m).name;
    results(m).params         = p_opt;
    results(m).h_wF           = p_opt(1:2);
    results(m).h_wR           = p_opt(3:4);
    results(m).padfrac_params = p_opt(5:end);
    results(m).rmse_F         = rmse;
    results(m).sse            = sse;
    results(m).nresid         = nResid;
    results(m).nparams        = kParams;
    results(m).aicc           = aicc;

    fprintf('  h_wF:  x1 = %.4f, b1 = %.4f\n', p_opt(1), p_opt(2));
    fprintf('  h_wR:  x1 = %.4f, b1 = %.4f\n', p_opt(3), p_opt(4));
    fprintf('  PadFrac params: %s\n', mat2str(p_opt(5:end), 5));
    fprintf('  RMSE = %.2f degF | AICc = %.1f\n', rmse, aicc);
end

%% ================== COMPARE MODELS ==================
fprintf('\n================ MODEL COMPARISON ================\n');
fprintf('%-38s %8s %8s %8s\n', 'Model', '#params', 'RMSE(F)', 'AICc');
for m = 1:numel(results)
    fprintf('%-38s %8d %8.2f %8.1f\n', results(m).name, results(m).nparams, results(m).rmse_F, results(m).aicc);
end
[~, best_idx] = min([results.aicc]);
fprintf('\n>>> Recommended model (lowest AICc, balances fit vs. complexity): %s <<<\n', results(best_idx).name);
fprintf('    If you just want the lowest raw error regardless of overfit risk, compare RMSE instead.\n');

%% ================== PLOT BEST MODEL: MEASURED VS SIMULATED ==================
best = results(best_idx);
best_fun = models(best_idx).fun;

for k = 1:nFiles
    ds = datasets(k);
    predF_front = run_sim_opt(ds.t, ds.velx, ds.frontpressure, ds.fr_temp_F, ds.Tbias_brake, ...
        best.h_wF(1), best.h_wF(2), best_fun, best.padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.fl_omega_wheel, ds.fr_omega_wheel, VehicleMass, RotorMass_front, RotorArea_front, I, WheelR, TambK);
    predF_rear = run_sim_opt(ds.t, ds.velx, ds.rearpressure, ds.rr_temp_F, 1 - ds.Tbias_brake, ...
        best.h_wR(1), best.h_wR(2), best_fun, best.padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.rl_omega_wheel, ds.rr_omega_wheel, VehicleMass, RotorMass_rear, RotorArea_rear, I, WheelR, TambK);

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
fprintf('x1f = %.6f;   %% h_wF slope\n', best.h_wF(1));
fprintf('b1f = %.6f;   %% h_wF intercept\n', best.h_wF(2));
fprintf('x1r = %.6f;   %% h_wR slope\n', best.h_wR(1));
fprintf('b1r = %.6f;   %% h_wR intercept\n', best.h_wR(2));
fprintf('PadFrac params (%s):\n', models(best_idx).name);
for i = 1:numel(best.padfrac_params)
    fprintf('  p(%d) = %.8g\n', i, best.padfrac_params(i));
end
fprintf('\nNOTE: if the best model has more than 2 PadFrac params (i.e. is not\n');
fprintf('the plain linear-in-T baseline), you must update the PadFrac formula\n');
fprintf('inside run_sim (in brake_temp_sim.m) to match models(%d).fun above -\n', best_idx);
fprintf('the old hardcoded "PadFrac = prevTemp*x2 + b2" line will not use the\n');
fprintf('new pressure-dependent terms.\n');


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
    I, WheelR, TambK)
% Pools front + rear residuals across ALL datasets into one vector for
% lsqnonlin. x1,b1 (h_w) and the PadFrac params are shared/global across
% every dataset and both corners - only Tbias flips between front/rear.

x1f_p = p(1); b1f_p = p(2); x1r_p = p(3); b1r_p = p(4); padfrac_params = p(5:end);

residual = [];
for k = 1:numel(datasets)
    ds = datasets(k);

    predF_front = run_sim_opt(ds.t, ds.velx, ds.frontpressure, ds.fr_temp_F, ds.Tbias_brake, ...
        x1f_p, b1f_p, padfrac_fun, padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.fl_omega_wheel, ds.fr_omega_wheel, VehicleMass, RotorMass_front, RotorArea_front, I, WheelR, TambK);

    predF_rear = run_sim_opt(ds.t, ds.velx, ds.rearpressure, ds.rr_temp_F, 1 - ds.Tbias_brake, ...
        x1r_p, b1r_p, padfrac_fun, padfrac_params, ds.total_regen_power, ds.Edrag, ...
        ds.rl_omega_wheel, ds.rr_omega_wheel, VehicleMass, RotorMass_rear, RotorArea_rear, I, WheelR, TambK);

    residual = [residual; predF_front(:) - ds.fr_temp_F(:); predF_rear(:) - ds.rr_temp_F(:)]; %#ok<AGROW>
end

residual(~isfinite(residual)) = 0;  % safety net against a rare divergent step
end


function RotorTempArrayF = run_sim_opt(t, velx, BrakePress, brakeTempArray, Tbias, ...
    x1_p, b1_p, padfrac_fun, padfrac_params, total_regen_power, Edrag, ...
    omega_wheel_L, omega_wheel_R, VehicleMass, RotorMass, RotorArea, I, WheelR, TambK) %#ok<INUSD>
% Same physics as run_sim() in brake_temp_sim.m, except PadFrac is now
% evaluated by an arbitrary function of (rotor temp [K], applied
% pressure [psi]) instead of a hardcoded linear-in-T formula.

RotorTempArrayK = zeros(size(t));
RotorTempArrayK(1) = (brakeTempArray(1) - 32) * (5/9) + 273.15;
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
