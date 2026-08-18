%% ================================================================
%  Brake Pad Specific Power (Heat Flux) & Fade-Onset Analysis
% ================================================================
%  Uses FIXED, already-calibrated rotor cooling and pad-fraction
%  coefficients (see "PADFRAC MODEL DEFINITIONS" below) to walk one or
%  more driving datasets and compute, for every identified braking event:
%
%     q'' * A_pad = (K_Linear + K_rot - D_aero) / t_decel   [applied per-pad]
%
%  where q'' is the average heat flux (W/cm^2) the pad experiences
%  during that event. PadFrac is applied to route the appropriate
%  share of that energy into the PAD (as opposed to the rotor), and
%  the same regen-energy subtraction, front/rear Tbias split, and
%  0.82 / 0.8 derating factors used in brake_temp_sim.m are applied
%  for full consistency with the calibrated thermal model.
%
%  Also empirically detects pad fade onset: a dynamics-based (mu-
%  independent) actual brake torque is computed per corner from wheel
%  deceleration, compared against a pressure+nominal-mu PREDICTED
%  torque, and the ratio (mu_actual/mu_nominal) is tracked against
%  specific power and temperature to find where the pads start
%  measurably underperforming their nominal friction coefficient.
%
%  Supports evaluating MULTIPLE padFrac models in one run (see
%  "PADFRAC MODEL DEFINITIONS") - every output below is generated once
%  per model, plus a final cross-model summary figure comparing specific
%  power results across models, to show how much the choice of padFrac
%  model actually impacts the results.
%
%  Outputs (generated once per padFrac model):
%    1) Specific power (W/cm^2) vs. time-since-event-start, for every
%       braking event found, front and rear pads plotted separately.
%    2) The average specific power across all events (front & rear).
%    3) A box-and-whisker plot of each event's average specific power
%       (front vs. rear).
%    4) Heatmaps of mu_actual/mu_nominal vs. specific power and
%       temperature, for both instantaneous and lagged flux.
%    5) The empirically-derived fade-onset specific power, front/rear.
%    6) Probability distributions of average specific power by drive
%       type (endurance, autocross, other/misc), front/rear overlaid,
%       with a table identifying each distribution's best-fit type and
%       summary statistics.
%    7) Scatter plots correlating specific power against deceleration
%       rate, starting speed, regen percentage, event duration, and
%       rotor temperature at event start.
%    8) A regen-braking-failure worst-case margin check: back-calculates
%       the pad area needed to stay under the fade-onset flux.
%  Plus, after all models are evaluated: a cross-model comparison figure
%  (box-and-whisker + overlaid distributions of specific power per model).
%
%  Datasets are selected via a multi-select file picker, same
%  21/22/27-column auto-detected formats as BrakeDataAnalysis.m.
% ================================================================

clc; clear; close all
cfg.SkipTimeCropPrompt = true;

%% ================== VEHICLE / BRAKE CONSTANTS ==================
% All vehicle/brake/sensor parameters come from the CENTRALIZED spreadsheet
% VehicleParameters.xlsx so this script and BrakeCoeffOptimizer.m can never
% drift apart. Every sizing parameter there is declared strictly PER SINGLE
% COMPONENT (one pad / one rotor / one corner / one caliper) - see
% loadVehicleParams.m's header for the full convention and for why the
% piston-area and pad-area scalings live in the loader.
VP = loadVehicleParams();

VehicleMass     = VP.VehicleMass;        % whole vehicle
RotorMass_front = VP.RotorMass_front;    % ONE rotor
RotorMass_rear  = VP.RotorMass_rear;     % ONE rotor
RotorArea_front = VP.RotorArea_front;    % ONE rotor
RotorArea_rear  = VP.RotorArea_rear;     % ONE rotor
I               = VP.I_corner;           % ONE corner
gear_ratio      = VP.gear_ratio;
TambC           = VP.TambC_fallback;     % fallback; per-dataset ambient derived at load time
BrakeFrac       = VP.BrakeFrac;
CalibrationFactor = VP.CalibrationFactor;

% Pad areas -- PER SINGLE PAD. A corner has `pads_per_corner` of these
% rubbing one rotor, so corner-level pad energy must be divided by
% pads_per_corner before dividing by this area.
A_pad_front_cm2 = VP.A_pad_front_cm2;
A_pad_rear_cm2  = VP.A_pad_rear_cm2;
pads_per_corner = VP.pads_per_corner;

velx_threshold = VP.velx_threshold;  % m/s, cleaning threshold
min_pressure   = VP.min_pressure;    % psi, threshold to consider a corner's brake engaged

% NOTE: the calibrated padFrac model(s) - cooling coefficients and
% PadFrac(T,P) formula(s) - are defined further below, in "PADFRAC MODEL
% DEFINITIONS", after the datasets are loaded (two of the default models
% need data-driven normalization constants computed from the loaded data).

%% ================== EVENT-DETECTION SETTINGS ==================
% A vehicle-level "braking event" = (frontpressure>min_pressure OR
% rearpressure>min_pressure) AND decelerating. Short gaps are bridged
% and very short/noisy blips are dropped -- see conversation notes:
%   - bridges brief pressure dips (trail-braking / modulation) so one
%     stop isn't split into fragments
%   - drops single-sample ADC noise spikes
%   - does NOT fix: pure-regen stops with pressure that never crosses
%     min_pressure (correctly excluded - no pad heating anyway)
bridge_gap_s     = 0.2;   % bridge gaps in active-braking shorter than this
min_event_dur_s  = 0.2;   % drop merged events shorter than this

%% ================== FADE-ONSET DETECTION SETTINGS ==================
% "Actual" torque comes from wheel/vehicle dynamics (KE loss, regen- and
% aero-corrected) - it does NOT depend on the mu_temp_table at all.
% "Predicted" torque comes from pressure * piston area * nominal mu (the
% SAME formula already used to build Tbias_brake). Their ratio is the
% fade signal: mu_ratio = T_actual / T_predicted = mu_actual / mu_nominal.
min_omega_wheel_rad_s = 3.0;   % below this, Power/omega -> torque is noisy; exclude
min_pressure_muratio_psi = 150; % psi - stricter than the general min_pressure engagement
                                % gate (5 psi). T_predicted ~ pressure * piston area *
                                % nominal mu * radius, so right at min_pressure it's only a
                                % few N*m - too small for the independently-computed
                                % T_actual's finite-difference noise to divide by cleanly.
                                % Mirrors min_omega_wheel_rad_s's role on the numerator side.
                                % CALIBRATED against curated_8-2 (14 files): a pressure-
                                % threshold sweep (5 to 1000 psi) showed the bulk of the
                                % mu_ratio distribution (p99) falls from ~19x at 5 psi to
                                % ~1.5x/2.0x (front/rear) at 150 psi, while retaining ~75-80%
                                % of samples (well above what analyzeFadeOnset's
                                % min_samples_per_bin needs). A handful of individual outlier
                                % points (up to ~7-13x) persist at ALL pressure levels tested,
                                % including high pressure - these are NOT a low-pressure/
                                % small-denominator artifact and this gate does not remove
                                % them; they most likely trace to single-timestep timing
                                % mismatch between the instantaneous pressure reading and the
                                % finite-difference wheel-deceleration signal. They are left
                                % visible (unclipped) in Figure 3 by design - see plotFadeHeatmap.
mu_ratio_fade_threshold = 0.9; % flag fade onset where mu_ratio drops below this
n_temp_bins      = 4;    % quantile bins used to control for temperature
n_flux_bins      = 7;    % quantile bins of specific power within each temp bin
min_samples_per_bin = 15; % ignore bins with too few samples to trust
lagged_flux_window_s = 1.0;  % trailing moving-average window for q_lag

%% ================== REGEN-FAILURE WORST-CASE ASSUMPTIONS ==================
% >>> These are placeholder engineering assumptions - review/edit before
%     trusting the back-calculated required pad area. <<<
% Worst case modeled as: vehicle at RegenFailure_Speed_mps, 100% of the
% stop handled by friction brakes (zero regen), decelerating at a
% constant RegenFailure_Decel_g, front/rear split per
% RegenFailure_TbiasFront (defaults to the mean Tbias_brake observed in
% the loaded data - edit if you want a specific/worst-case bias instead).
% Peak specific power for a constant-deceleration stop occurs at t=0
% (highest speed), so that's what's compared against the fade threshold.
RegenFailure_Speed_mps    = [];    % [] = use the max speed observed in loaded data
RegenFailure_Decel_g      = 1.5;   % assumed achievable peak deceleration
RegenFailure_TbiasFront   = [];    % [] = use the mean observed Tbias_brake

%% ================== SELECT DATASETS ==================
[files, path] = uigetfile('*.txt', 'Select driving dataset file(s)', 'MultiSelect', 'on');
if isequal(files, 0)
    error('No files selected.');
end
if ischar(files)
    files = {files};
end
nFiles = numel(files);

% The loaded parameter struct P is passed straight through to
% compute_derived_quantities - no intermediate aeroParams/gearParams copies,
% so there is exactly one source for every constant.

TambK = TambC + 273.15;   % fallback only; each dataset derives its own below

%% ================== LOAD & CACHE PER-FILE DATA (model-independent) ==================
% File I/O, parsing, derived-quantity computation, drive-type detection,
% and braking-event identification all depend only on the raw telemetry
% and vehicle/event-detection constants above - NOT on which padFrac
% model is being evaluated. So this all runs ONCE here, cached into
% `datasets`, regardless of how many padFrac models are evaluated below
% (only the simulation itself, which does depend on padFrac, is repeated
% per model in the loop further down).
datasets = struct('name', {}, 't', {}, 'derived', {}, 'Edrag', {}, 'driveType', {}, ...
    'event_ranges', {}, 'TambK', {});

maxSpeedSeen = 0;
sumTbias = 0; countTbias = 0;

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
    derived = compute_derived_quantities(parsed, VP);
    Edrag   = compute_edrag(parsed.t, derived.velx, derived.F_aero);
    driveType = detectDriveType(fname);
    fprintf('  Drive type: %s\n', driveType);

    % Per-dataset ambient, derived the SAME way BrakeCoeffOptimizer.m does
    % when it fits the h_w/PadFrac coefficients: average of this dataset's
    % first front and rear rotor readings. Previously this script applied
    % those fitted coefficients against a hardcoded 25 C, i.e. outside the
    % conditions they were fit in.
    firstFrontF = derived.fr_temp_F(1);
    firstRearF  = derived.rr_temp_F(1);
    if isfinite(firstFrontF) && isfinite(firstRearF)
        TambC_dataset = ((firstFrontF + firstRearF) / 2 - 32) * (5/9);
    elseif isfinite(firstFrontF)
        TambC_dataset = (firstFrontF - 32) * (5/9);
    elseif isfinite(firstRearF)
        TambC_dataset = (firstRearF - 32) * (5/9);
    else
        TambC_dataset = TambC;
        warning('  Dataset %s: no valid first rotor temps - falling back to %.1f C.', files{k}, TambC);
    end
    fprintf('  Dataset ambient: %.1f C\n', TambC_dataset);

    maxSpeedSeen = max(maxSpeedSeen, max(derived.velx, [], 'omitnan'));
    validTbias = derived.Tbias_brake(isfinite(derived.Tbias_brake) & derived.Tbias_brake > 0);
    sumTbias = sumTbias + sum(validTbias);
    countTbias = countTbias + numel(validTbias);

    % Vehicle-level braking-event ranges - identical across every padFrac
    % model (depends only on pressure/velocity, not PadFrac), so computed
    % once here rather than redundantly once per model.
    DS = [0; diff(derived.velx)];
    DS(DS < -2.5) = 0;  % matches the >2.5g exclusion used in the sim
    active = (derived.frontpressure > min_pressure | derived.rearpressure > min_pressure) & (DS < 0);
    event_ranges = find_events(active, parsed.t, bridge_gap_s, min_event_dur_s);
    fprintf('  Found %d braking event(s).\n', size(event_ranges, 1));

    datasets(end+1) = struct('name', files{k}, 't', parsed.t, 'derived', derived, ...
        'Edrag', Edrag, 'driveType', driveType, 'event_ranges', event_ranges, ...
        'TambK', TambC_dataset + 273.15); %#ok<SAGROW>
end

if isempty(datasets)
    error('No usable datasets loaded.');
end

%% ================== DATA-DRIVEN PADFRAC ANCHOR CONSTANTS (models 5 & 6) ==================
% Two of the six default padFrac models below (5: "Anchored to
% effusivity-based ideal", 6: "Saturating logistic in T,P") were fit in
% BrakeCoeffOptimizer.m against normalization constants (Tmid_K, dT,
% Pmid, dP) derived from that run's loaded dataset's observed T/P range,
% plus a fixed material-property anchor (PadFrac_ideal). None of these
% were logged as standalone numbers alongside the fitted coefficients, so
% they're recomputed here, replicating BrakeCoeffOptimizer.m's exact
% formulas, from the SAME curated_8-2 files those coefficients were fit
% against. If this script is later pointed at a different dataset,
% models 5/6 will drift from their logged fit (models 1-4 are unaffected,
% they don't reference these constants).
allT_K = []; allP = [];
for k = 1:numel(datasets)
    d = datasets(k).derived;
    allT_K = [allT_K; (d.fr_temp_F - 32)*(5/9) + 273.15; (d.rr_temp_F - 32)*(5/9) + 273.15]; %#ok<AGROW>
    allP   = [allP; d.frontpressure; d.rearpressure]; %#ok<AGROW>
end
Tmin_K = min(allT_K); Tmax_K = max(allT_K);
dT = max(Tmax_K - Tmin_K, 1);
Pmax = max(allP);
dP = max(Pmax, 1);
Tmid_K = (Tmin_K + Tmax_K) / 2;
Pmid   = Pmax / 2;

k_rotor = 42.7; rho_rotor = 7850; cp_rotor = 477;       % 4130 steel (BrakeCoeffOptimizer.m)
xi_rotor = sqrt(k_rotor * rho_rotor * cp_rotor);
k_pad_mid = 3.0; rho_pad_mid = 2200; cp_pad_mid = 1000; % Porterfield R4-1 estimate (BrakeCoeffOptimizer.m)
xi_pad_mid = sqrt(k_pad_mid * rho_pad_mid * cp_pad_mid);
PadFrac_ideal = xi_pad_mid / (xi_pad_mid + xi_rotor);

%% ================== PADFRAC MODEL DEFINITIONS ==================
% The 6 candidate padFrac models, REFIT against carData\curated_8-2 (14
% files) on 14-Aug-2026 under the corrected per-component physics (see
% SIZING_AUDIT.md). The previous coefficients - fit before the rotational-KE
% aggregation, per-corner regen attribution, and per-dataset ambient fixes -
% no longer correspond to this model and were replaced wholesale. All values
% are the 'joint' fitting strategy, which beat 'twostage' on AICc for every
% model. Each model carries its OWN cooling coefficients (x1f/b1f/x1r/b1r)
% since these were jointly fit per padFrac model - reusing one shared
% cooling fit across all 6 would confound the "how much does padFrac model
% choice matter" comparison this array exists for.
%
% NOTE: coefficients below carry ~5 significant figures (the precision
% BrakeCoeffOptimizer prints per model). Re-run it if you need more.
padFracModels(1).name = 'Linear in T only';
padFracModels(1).x1f = 2.366902; padFracModels(1).b1f = 21.133072;
padFracModels(1).x1r = 3.641386; padFracModels(1).b1r = 29.359405;
padFracModels(1).fun = @(T,P) 0.00046090397.*T + 0.13389237;

% padFracModels(2).name = 'Linear, independent T and P';
% padFracModels(2).x1f = 2.3669; padFracModels(2).b1f = 21.1331;
% padFracModels(2).x1r = 3.6414; padFracModels(2).b1r = 29.3594;
% padFracModels(2).fun = @(T,P) 0.00046091.*T + 2.3374e-14.*P + 0.13389;
% 
% padFracModels(3).name = 'Linear with T*P interaction';
% padFracModels(3).x1f = 0.2886; padFracModels(3).b1f = 23.6252;
% padFracModels(3).x1r = 1.7692; padFracModels(3).b1r = 30.9199;
% padFracModels(3).fun = @(T,P) 0.001538.*T + 7.2797e-13.*P + (-9.3245e-07).*T.*P + (-0.057921);
% 
% padFracModels(4).name = 'Quadratic in T, linear in P';
% padFracModels(4).x1f = 2.4898; padFracModels(4).b1f = 21.4565;
% padFracModels(4).x1r = 3.9885; padFracModels(4).b1r = 26.7070;
% padFracModels(4).fun = @(T,P) -0.0022399.*T + 2.2497e-06.*T.^2 + 2.3373e-14.*P + 0.8988;
% 
% padFracModels(5).name = 'Anchored to effusivity-based ideal';
% padFracModels(5).x1f = 3.7213; padFracModels(5).b1f = 22.6092;
% padFracModels(5).x1r = 4.9567; padFracModels(5).b1r = 30.1357;
% padFracModels(5).fun = @(T,P) PadFrac_ideal + 0.13218.*(T-Tmid_K)/dT + (-0.15).*(P-Pmid)/dP + 0.05;
% 
% padFracModels(6).name = 'Saturating logistic in T,P';
% padFracModels(6).x1f = 1.7024; padFracModels(6).b1f = 23.0513;
% padFracModels(6).x1r = 2.6705; padFracModels(6).b1r = 30.5203;
% padFracModels(6).fun = @(T,P) 0.58619 ./ (1 + exp(-(4.1327.*(T-Tmid_K)/dT + (-5.2435).*(P-Pmid)/dP + 0.45673)));

C = struct('VehicleMass', VehicleMass, 'RotorMass_front', RotorMass_front, 'RotorMass_rear', RotorMass_rear, ...
    'RotorArea_front', RotorArea_front, 'RotorArea_rear', RotorArea_rear, 'I', I, ...
    'A_pad_front_cm2', A_pad_front_cm2, 'A_pad_rear_cm2', A_pad_rear_cm2, ...
    'pads_per_corner', pads_per_corner, ...
    'BrakeFrac', BrakeFrac, 'CalibrationFactor', CalibrationFactor, 'min_pressure', min_pressure, ...
    'min_omega_wheel_rad_s', min_omega_wheel_rad_s, 'min_pressure_muratio_psi', min_pressure_muratio_psi, ...
    'mu_ratio_fade_threshold', mu_ratio_fade_threshold, 'n_temp_bins', n_temp_bins, 'n_flux_bins', n_flux_bins, ...
    'min_samples_per_bin', min_samples_per_bin, 'lagged_flux_window_s', lagged_flux_window_s, ...
    'RegenFailure_Speed_mps', RegenFailure_Speed_mps, 'RegenFailure_Decel_g', RegenFailure_Decel_g, ...
    'RegenFailure_TbiasFront', RegenFailure_TbiasFront, ...
    'maxSpeedSeen', maxSpeedSeen, 'sumTbias', sumTbias, 'countTbias', countTbias);

%% ================== RUN ANALYSIS PER PADFRAC MODEL ==================
modelResults = struct('name', {}, 'avg_q_front', {}, 'avg_q_rear', {});
for m = 1:numel(padFracModels)
    modelResults(m) = run_model_analysis(datasets, padFracModels(m), m, C);
end

%% ================== CROSS-MODEL SUMMARY: HOW MUCH DOES PADFRAC MODEL CHOICE MATTER? ==================
fprintf('\n================ CROSS-MODEL SUMMARY ================\n');
for m = 1:numel(modelResults)
    fprintf('M%d: %s\n', m, modelResults(m).name);
end

figure('Name', 'Cross-Model Specific Power Comparison');

% --- Left: grouped box-and-whisker, one box per model x corner ---
subplot(1,2,1);
allVals = []; grp = {}; order = {};
for m = 1:numel(modelResults)
    fLabel = sprintf('M%d Front', m); rLabel = sprintf('M%d Rear', m);
    order = [order, {fLabel, rLabel}]; %#ok<AGROW>
    allVals = [allVals; modelResults(m).avg_q_front(:); modelResults(m).avg_q_rear(:)]; %#ok<AGROW>
    grp = [grp; repmat({fLabel}, numel(modelResults(m).avg_q_front), 1); ...
                 repmat({rLabel}, numel(modelResults(m).avg_q_rear), 1)]; %#ok<AGROW>
end
boxplot(allVals, grp, 'GroupOrder', order);
ylabel('Avg Specific Power per Event (W/cm^2)');
title('Spread by Model (box = IQR, line = median)');
xtickangle(45);
grid on;

% --- Right: overlaid fitted PDFs, one per model (front+rear pooled) ---
subplot(1,2,2); hold on; grid on;
colors = lines(numel(modelResults));
for m = 1:numel(modelResults)
    pooled = [modelResults(m).avg_q_front(:); modelResults(m).avg_q_rear(:)];
    distFit = fit_distribution_summary(pooled);
    if ~isempty(distFit.pd)
        xg = linspace(min(pooled), max(pooled), 200);
        plot(xg, pdf(distFit.pd, xg), 'Color', colors(m,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('M%d: %s', m, modelResults(m).name));
    end
end
xlabel('Avg Specific Power (W/cm^2), front+rear pooled'); ylabel('Probability Density');
title('Distribution Shift by Model');
legend('Location', 'best', 'FontSize', 7, 'Interpreter', 'none');


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


function derived = compute_derived_quantities(parsed, VP)
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

% Clamp force uses the SINGLE-SIDE piston area (loadVehicleParams halves the
% per-caliper piston count for opposed calipers). The caliper is a closed
% force loop: the far side reacts the load, it does not add to it. The
% factor of 2 in the torque expression below is the two friction FACES, not
% the two sides' pistons - applying both would double-count.
front_piston_area_side = VP.front_piston_area_per_side;
rear_piston_area_side  = VP.rear_piston_area_per_side;

% Friction lever arm is the MEAN PAD RADIUS, an explicit documented input -
% not the rotor outer radius, which overstates the moment arm.
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

F_aero_open   = VP.aero_open_a   * speed_mph.^2 + VP.aero_open_b   * speed_mph + VP.aero_open_c;
F_aero_closed = VP.aero_closed_a * speed_mph.^2 + VP.aero_closed_b * speed_mph + VP.aero_closed_c;
F_aero = F_aero_open .* double(parsed.drs_state == 1) + F_aero_closed .* double(parsed.drs_state == 0);
F_aero = max(F_aero, 0);

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
derived.fl_omega_wheel    = fl_omega_wheel;
derived.fr_omega_wheel    = fr_omega_wheel;
derived.rl_omega_wheel    = rl_omega_wheel;
derived.rr_omega_wheel    = rr_omega_wheel;
% Regen is measured PER CORNER, so it is kept per axle rather than pooled to
% a single vehicle-level total. Pooling it (the previous behavior) and
% subtracting before the Tbias/corner split let front regen offset REAR pad
% heat in proportion to brake bias, which is not physical.
derived.front_axle_regen_power = fl_regen_power + fr_regen_power;
derived.rear_axle_regen_power  = rl_regen_power + rr_regen_power;
derived.total_regen_power      = derived.front_axle_regen_power + derived.rear_axle_regen_power;
% Per-corner PREDICTED (single-side clamp force x nominal mu x mean pad
% radius x 2 friction faces) torque - exposed for the dynamics-vs-predicted
% fade comparison. Sign is dropped (magnitude only) since we only need it
% for a ratio against the (also-positive) actual torque. Front/rear each
% have a single hydraulic circuit, so fl==fr and rl==rr by construction.
derived.T_predicted_front = abs(fl_Tbrake);
derived.T_predicted_rear  = abs(rl_Tbrake);
derived.mu_nominal_front  = mu_front;
derived.mu_nominal_rear   = mu_rear;
end


function Edrag = compute_edrag(t, velx, F_aero)
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


function sim = simulate_pad_power(t, velx, BrakePress, Tbias, x1_p, b1_p, padfrac_fun, ...
    axle_regen_power, Edrag, omega_wheel_L, omega_wheel_R, omega_all, ...
    VehicleMass, RotorMass, RotorArea, I, TambK, A_pad_cm2, pads_per_corner, BrakeFrac, CalibrationFactor, ...
    min_pressure, corner_split, T_rotor_start_K, T_predicted, min_omega_wheel_rad_s, min_pressure_muratio)
% Walks the full timeseries computing rotor temp (needed for PadFrac's
% temperature dependence), per-step pad energy, and instantaneous
% specific power (W/cm^2). PadFrac uses the SAME (1-PadFrac)/PadFrac
% split convention as brake_temp_sim.m: PadFrac is the pad's own share,
% (1-PadFrac) is the rotor's share.
%
% BASIS: every sizing input is per single component. The energy chain is
%   whole vehicle -> x Tbias -> one axle -> x corner_split -> ONE CORNER
% and a corner contains ONE rotor but `pads_per_corner` PADS. So the rotor
% branch divides corner energy by one RotorMass directly, while the pad
% branch must additionally divide by pads_per_corner before dividing by
% A_pad_cm2 (which is ONE pad's area). Omitting that step is what made
% specific power read 2x high previously.
%
% Also computes, per step, a DYNAMICS-BASED actual brake torque (mu-
% independent - derived from wheel/vehicle KE loss, not from pressure and
% an assumed mu) and the resulting mu_ratio = T_actual/T_predicted fade
% signal. Unlike pad_energy/q_inst (which are deliberately discounted by
% BrakeFrac/CalibrationFactor/PadFrac - empirical THERMAL calibration
% factors), T_actual is NOT discounted by those factors: it is a
% mechanical torque estimate, and those three factors describe how
% friction energy partitions into heat afterward, not how much torque
% the brake actually applied. Aero and regen ARE still subtracted, since
% both physically bypass the friction brake entirely.
%
% `omega_all` is a 4-column [fl fr rl rr] matrix of per-corner wheel speeds,
% used to sum rotational KE strictly per corner rather than scaling one
% axle's mean omega to all four wheels.
%
% `axle_regen_power` is THIS axle's regen only (not the vehicle total), so
% regen is attributed to the corners that actually produced it.
%
% Also exposes per-step regen_energy/friction_energy at the CORNER level
% (both zero outside the active-braking branch, same convention as
% pad_energy) so callers can compute regen's share of an event's braking
% energy without recomputing this loop's KE/aero logic.

n = length(t);
RotorTempArrayK = zeros(n, 1);
RotorTempArrayK(1) = T_rotor_start_K;   % measured first rotor temp (see caller)
pad_energy      = zeros(n, 1);   % J, per SINGLE PAD, per step
q_inst          = zeros(n, 1);   % W/cm^2, per SINGLE PAD, per step
T_actual        = nan(n, 1);     % Nm, per single corner, dynamics-based
mu_ratio        = nan(n, 1);     % T_actual / T_predicted
regen_energy    = zeros(n, 1);   % J, per step, THIS CORNER
friction_energy = zeros(n, 1);   % J, per step, THIS CORNER (post aero/regen)

for i = 2:n
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
    PadFrac  = max(min(padfrac_fun(prevTemp, BrakePress(i)), 1), 0);
    SpecHeat = (0.0005 * prevTemp + 0.2813) * 1000;

    if DS < -2.5
        DS = 0;
    end

    if DS < 0 && BrakePress(i) > min_pressure
        Energy1 = 0.5 * VehicleMass * (prevSpeed^2 - newSpeed^2);
        % Rotational KE summed strictly PER CORNER: I is one corner's
        % inertia, so each corner contributes 0.5*I*(w_prev^2 - w_new^2)
        % using its OWN omega. (Previously this squared one axle's MEAN
        % omega and scaled it x4, which both mixes axles and relies on
        % mean(w)^2 == mean(w^2).)
        Energy2 = sum(0.5 * I * (omega_all(i-1,:).^2 - omega_all(i,:).^2));
        Energy  = Energy1 + Energy2;

        omegaN  = (omega_wheel_L(i) + omega_wheel_R(i)) / 2;   % this axle, for torque
        AeroFrac  = min(Edrag(i) / max(Energy, 1), 1);

        % Split the vehicle's KE loss down to THIS CORNER first, then remove
        % aero and this corner's own regen. Regen used to be pooled across
        % all four corners and subtracted before the split, which let front
        % regen offset rear pad heat in proportion to brake bias.
        Energy_corner    = Energy * Tbias(i) * corner_split * (1 - AeroFrac);
        regen_energy(i)  = abs(axle_regen_power(i)) * tbrake * corner_split;
        friction_energy(i) = max(Energy_corner - regen_energy(i), 0);

        % Rotor share - one corner has ONE rotor, so no further division
        CorrectedEnergyRotor = friction_energy(i) * (1 - PadFrac) * BrakeFrac * CalibrationFactor;
        deltaTK = CorrectedEnergyRotor / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = deltaTK + prevTemp;

        qout = (RotorTempArrayK(i) - TambK) / Rrotor;
        Eout = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = RotorTempArrayK(i) - deltaTKout;

        % Pad share - one corner has `pads_per_corner` pads sharing this
        % heat, and A_pad_cm2 is ONE pad's area, so divide by the count.
        CorrectedEnergyPad = friction_energy(i) * PadFrac * BrakeFrac * CalibrationFactor;
        pad_energy(i) = CorrectedEnergyPad / pads_per_corner;
        q_inst(i)     = pad_energy(i) / tbrake / A_pad_cm2;

        % Dynamics-based (mu-independent) actual torque: mechanical brake
        % power at this corner, divided by wheel angular velocity. Gated
        % on a minimum omega to avoid Power/omega blowing up near a stop.
        if omegaN > min_omega_wheel_rad_s
            T_actual(i) = (friction_energy(i) / tbrake) / omegaN;
            if isfinite(T_predicted(i)) && T_predicted(i) > 1e-6 && BrakePress(i) > min_pressure_muratio
                mu_ratio(i) = T_actual(i) / T_predicted(i);
            end
        end
    else
        qout = (prevTemp - TambK) / Rrotor;
        Eout = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = prevTemp - deltaTKout;
        % pad_energy(i), q_inst(i), T_actual(i), mu_ratio(i), regen_energy(i),
        % friction_energy(i) remain 0/NaN
    end
end

sim.RotorTempArrayK  = RotorTempArrayK;
sim.pad_energy       = pad_energy;
sim.q_inst           = q_inst;
sim.T_actual         = T_actual;
sim.mu_ratio         = mu_ratio;
sim.regen_energy     = regen_energy;
sim.friction_energy  = friction_energy;
end


function event_ranges = find_events(active, t, bridge_gap_s, min_event_dur_s)
% Run-length-encodes the `active` logical vector into [start_idx end_idx]
% ranges, bridges gaps shorter than bridge_gap_s, then drops events
% shorter than min_event_dur_s.

n = length(active);
raw = [];  % [start end]
i = 1;
while i <= n
    if active(i)
        j = i;
        while j < n && active(j+1)
            j = j + 1;
        end
        raw = [raw; i, j]; %#ok<AGROW>
        i = j + 1;
    else
        i = i + 1;
    end
end

if isempty(raw)
    event_ranges = zeros(0, 2);
    return
end

% Bridge gaps
merged = raw(1, :);
for k = 2:size(raw, 1)
    gap_time = t(raw(k,1)) - t(merged(end,2));
    if gap_time <= bridge_gap_s + 1e-9   % epsilon guards against fp rounding at the boundary
        merged(end, 2) = raw(k, 2);   % extend current event
    else
        merged = [merged; raw(k, :)]; %#ok<AGROW>
    end
end

% Drop too-short events
durations = t(merged(:,2)) - t(merged(:,1));
event_ranges = merged(durations >= min_event_dur_s, :);
end


function driveType = detectDriveType(fname)
% Primary signal: the FILENAME itself. Checked first because
% RawFormatOnly exports from F34BrakeDataPreprocessor.m always write
% "Segment type: raw_format_only" in their header regardless of the
% file's true drive type - the header can't distinguish an endurance
% rawformat file from an autocross one, but the filename usually can.
[~, baseName, ~] = fileparts(fname);
lowerName = lower(baseName);
if ~isempty(strfind(lowerName, 'endurance'))
    driveType = 'Endurance';
    return
end
if ~isempty(strfind(lowerName, 'autocross')) || ~isempty(strfind(lowerName, 'autox'))
    driveType = 'Autocross';
    return
end

% Fallback: the "# Segment type: <kind>" header, for files whose name
% doesn't happen to contain a recognizable keyword but WERE preprocessed
% normally (not RawFormatOnly).
driveType = '';
fid = fopen(fname, 'r');
if fid >= 0
    for lineNum = 1:10
        line = fgetl(fid);
        if ~ischar(line)
            break
        end
        tok = regexp(line, '#\s*Segment type:\s*(\S+)', 'tokens', 'once');
        if ~isempty(tok)
            kind = lower(tok{1});
            if ~isempty(strfind(kind, 'endurance'))
                driveType = 'Endurance';
            elseif ~isempty(strfind(kind, 'autocross'))
                driveType = 'Autocross';
            else
                driveType = 'Other/Misc';
            end
            break
        end
    end
    fclose(fid);
end

if isempty(driveType)
    options = {'Endurance', 'Autocross', 'Other/Misc'};
    [selIdx, ok] = listdlg('ListString', options, 'SelectionMode', 'single', ...
        'InitialValue', 3, 'Name', 'Drive Type', ...
        'PromptString', {'No "Segment type" header found in:', fname, ...
        'Select this file''s drive type:'});
    if ok
        driveType = options{selIdx};
    else
        driveType = 'Other/Misc';
    end
end
end



function plotFadeHeatmap(flux, tempF, muRatio, titleStr)
% 3D surface + scatter3, matching the visual style of BrakeCoeffOptimizer.m's
% "PadFrac vs Temperature(F) and Pressure" plot (surf + scatter3 + colorbar +
% fixed view angle) applied to mu_ratio vs (flux, temperature).
%
% Upper-tail outliers (Tukey IQR "extreme outlier" fence: > Q3 + 3*(Q3-Q1),
% computed from THIS panel's own mu_ratio values) are omitted from the plot
% so a handful of extreme single-sample spikes don't dominate the surface/
% color scale. The 3x (not the more common 1.5x "mild outlier") multiplier
% is deliberate: mu_ratio's distribution is heavily right-skewed with most
% mass compressed between 0 and ~1, so the 1.5x fence flags 6-9% of ALL
% samples - including the entire physically-meaningful 0.85-1.2 range right
% around the fade threshold - not just genuine spikes. 3x targets only the
% true outliers. Only the upper tail is ever filtered - low mu_ratio is
% genuine fade signal, never noise, so it's never excluded. This filtering
% is local to the plot only; it does NOT touch the fade_samples_*/
% analyzeFadeOnset data used by Output 5/8. The omitted count is reported
% both to the console and in the subplot title so the omission itself
% stays visible.
valid = isfinite(flux) & isfinite(tempF) & isfinite(muRatio);
flux = flux(valid); tempF = tempF(valid); muRatio = muRatio(valid);
nTotal = numel(muRatio);

if nTotal >= 4
    q = quantile(muRatio, [0.25, 0.75]);
    upperFence = q(2) + 3 * (q(2) - q(1));
    isOutlier = muRatio > upperFence;
else
    isOutlier = false(size(muRatio));
    upperFence = Inf;
end
nOutliers = sum(isOutlier);
fprintf('  %s: omitted %d/%d outlier point(s) (mu_ratio > %.3f)\n', titleStr, nOutliers, nTotal, upperFence);

flux = flux(~isOutlier); tempF = tempF(~isOutlier); muRatio = muRatio(~isOutlier);
titleStr = sprintf('%s (%d/%d omitted)', titleStr, nOutliers, nTotal);

if numel(flux) < 10 || range(flux) <= 0 || range(tempF) <= 0
    title([titleStr ' (insufficient data)']);
    axis off;
    return
end

flux_grid = linspace(min(flux), max(flux), 50);
temp_grid = linspace(min(tempF), max(tempF), 50);
[F_grid, T_grid] = meshgrid(flux_grid, temp_grid);

Z_grid = griddata(flux, tempF, muRatio, F_grid, T_grid, 'linear');
Z_grid = imgaussfilt(Z_grid, 1.5);

surf(F_grid, T_grid, Z_grid, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on;
scatter3(flux, tempF, muRatio, 18, muRatio, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.85);

colorbar;
xlabel('Specific Power (W/cm^2)');
ylabel('Rotor Temp (deg F, pad-face proxy)');
zlabel('mu_{actual} / mu_{nominal}');
title(titleStr, 'Interpreter', 'none');
view(45, 25);
grid on;
end


function result = analyzeFadeOnset(samples, threshold, n_temp_bins, n_flux_bins, min_samples_per_bin)
% samples: [q_inst, q_lag, T_F, mu_ratio] pooled per-sample data.
% Bins by temperature (quantiles) first to control for mu's own normal
% temperature dependence, then within each temperature bin bins by
% specific power (quantiles) and finds the first flux bin where mean
% mu_ratio drops below `threshold` and STAYS below it as flux keeps
% increasing (guards against flagging a single noisy dip). The overall
% onset is the median across temperature bins; the coefficient of
% variation (CV) across temperature bins is used elsewhere to judge
% which flux definition (instantaneous vs. lagged) correlates more
% consistently with the mu drop.

fluxCols  = [1, 2];
fluxNames = {'q_inst', 'q_lag'};

for fc = 1:numel(fluxCols)
    col = fluxCols(fc);
    flux   = samples(:, col);
    tempF  = samples(:, 3);
    muR    = samples(:, 4);

    valid = isfinite(flux) & isfinite(tempF) & isfinite(muR);
    flux = flux(valid); tempF = tempF(valid); muR = muR(valid);

    r = struct('onsetPerTempBin', [], 'onsetOverall', NaN, 'cv', NaN, 'nSamples', numel(flux));

    if numel(flux) >= min_samples_per_bin * 2 && range(tempF) > 0
        tempEdges = quantile(tempF, linspace(0, 1, n_temp_bins+1));
        tempEdges(1) = -inf; tempEdges(end) = inf;

        onsetPerBin = nan(n_temp_bins, 1);
        for tb = 1:n_temp_bins
            inTempBin = tempF > tempEdges(tb) & tempF <= tempEdges(tb+1);
            if sum(inTempBin) < min_samples_per_bin
                continue
            end
            fluxBin = flux(inTempBin);
            muBin   = muR(inTempBin);

            fluxEdges = unique(quantile(fluxBin, linspace(0, 1, n_flux_bins+1)));
            if numel(fluxEdges) < 3
                continue
            end
            nBinsActual = numel(fluxEdges) - 1;
            meanMuPerFluxBin  = nan(nBinsActual, 1);
            fluxBinCenters    = nan(nBinsActual, 1);

            for qb = 1:nBinsActual
                inFluxBin = fluxBin > fluxEdges(qb) & fluxBin <= fluxEdges(qb+1);
                if qb == 1
                    inFluxBin = inFluxBin | fluxBin == fluxEdges(qb);
                end
                if sum(inFluxBin) >= max(3, round(min_samples_per_bin/n_flux_bins))
                    meanMuPerFluxBin(qb) = mean(muBin(inFluxBin));
                    fluxBinCenters(qb)   = mean(fluxBin(inFluxBin));
                end
            end

            belowIdx = find(meanMuPerFluxBin < threshold);
            for bi = 1:numel(belowIdx)
                startBin = belowIdx(bi);
                remaining = meanMuPerFluxBin(startBin:end);
                remaining = remaining(isfinite(remaining));
                if ~isempty(remaining) && all(remaining < threshold)
                    onsetPerBin(tb) = fluxBinCenters(startBin);
                    break
                end
            end
        end

        r.onsetPerTempBin = onsetPerBin;
        validOnsets = onsetPerBin(isfinite(onsetPerBin));
        if ~isempty(validOnsets)
            r.onsetOverall = median(validOnsets);
            if r.onsetOverall ~= 0
                r.cv = std(validOnsets) / abs(mean(validOnsets));
            end
        end
    end

    result.(fluxNames{fc}) = r;
end
end


function report_fade_corner(label, fadeResult, A_cm2, avgQ, peakQ, lagWindow_s)
fprintf('--- %s ---\n', label);
fprintf('  Pad area: %.1f cm^2\n', A_cm2);
fprintf('  Current average specific power (all events): %.3f W/cm^2\n', avgQ);
fprintf('  Current peak per-event average specific power: %.3f W/cm^2\n', peakQ);

fluxNames = {'q_inst', 'q_lag'};
labels    = {'Instantaneous flux', sprintf('%.1fs-lagged flux', lagWindow_s)};
for i = 1:2
    r = fadeResult.(fluxNames{i});
    if isnan(r.onsetOverall)
        fprintf('  %s: insufficient data to determine a fade-onset threshold (n=%d samples).\n', ...
            labels{i}, r.nSamples);
    else
        fprintf('  %s fade-onset: %.3f W/cm^2 (n=%d samples, CV across temp bins = %.2f)\n', ...
            labels{i}, r.onsetOverall, r.nSamples, r.cv);
    end
end

if isfinite(fadeResult.q_inst.cv) && isfinite(fadeResult.q_lag.cv)
    if fadeResult.q_inst.cv <= fadeResult.q_lag.cv
        fprintf('  -> Instantaneous flux correlates more consistently with the mu drop across temperature bins.\n');
    else
        fprintf('  -> %s correlates more consistently with the mu drop across temperature bins.\n', labels{2});
    end
end
fprintf('\n');
end


function fitResult = fit_distribution_summary(data)
% Fits candidate distributions to `data` and selects the best by AIC
% (2*k - 2*loglikelihood), computed manually since fitdist objects don't
% expose a uniform AIC property across distribution types. Requires the
% Statistics and Machine Learning Toolbox (already an implicit dependency
% of this script via quantile/boxplot).
data = data(:);
data = data(isfinite(data));
fitResult = struct('n', numel(data), 'mean', NaN, 'median', NaN, 'std', NaN, ...
    'skew', NaN, 'bestName', 'insufficient data', 'bestParams', [], 'pd', []);

if fitResult.n == 0
    return
end

fitResult.mean   = mean(data);
fitResult.median = median(data);
fitResult.std    = std(data);
fitResult.skew   = skewness(data);

if fitResult.n < 5
    fitResult.bestName = sprintf('insufficient data (n=%d)', fitResult.n);
    return
end

candidates = {'Normal', 'Lognormal', 'Gamma', 'Weibull'};
bestAIC = Inf;
for c = 1:numel(candidates)
    try
        pd = fitdist(data, candidates{c});
        % Clip pdf values away from exact 0 before taking log: a single
        % far-outlier point can otherwise underflow pdf() to 0.0 and send
        % the log-likelihood to -Inf, disqualifying an otherwise-reasonable
        % fit via the isfinite(aic) check below instead of just penalizing
        % it (which is what should happen for one bad point).
        ll = sum(log(max(pdf(pd, data), realmin)));
        k  = numel(pd.ParameterNames);
        aic = 2*k - 2*ll;
        if isfinite(aic) && aic < bestAIC
            bestAIC = aic;
            fitResult.bestName   = candidates{c};
            fitResult.bestParams = pd.ParameterValues;
            fitResult.pd         = pd;
        end
    catch
        continue   % candidate not valid for this data (e.g. non-positive values)
    end
end
end


function [fitFront, fitRear] = plot_distribution_overlay(dataFront, dataRear, titleStr)
% Overlaid front (red) / rear (blue) probability-density histograms with
% each series' best-fit PDF curve drawn on top. Plots into the CURRENT
% axes (call subplot(...) beforehand), matching plotFadeHeatmap's style.
% Returns both fit results so the caller can also build a stats table
% from them without recomputing.
fitFront = fit_distribution_summary(dataFront);
fitRear  = fit_distribution_summary(dataRear);

hold on; grid on;
if fitFront.n > 0
    histogram(dataFront, 'Normalization', 'pdf', 'FaceColor', 'r', 'FaceAlpha', 0.35, ...
        'EdgeColor', 'none', 'DisplayName', 'Front');
end
if fitRear.n > 0
    histogram(dataRear, 'Normalization', 'pdf', 'FaceColor', 'b', 'FaceAlpha', 0.35, ...
        'EdgeColor', 'none', 'DisplayName', 'Rear');
end

if ~isempty(fitFront.pd)
    xg = linspace(min(dataFront), max(dataFront), 200);
    plot(xg, pdf(fitFront.pd, xg), 'r-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Front fit: %s', fitFront.bestName));
end
if ~isempty(fitRear.pd)
    xg = linspace(min(dataRear), max(dataRear), 200);
    plot(xg, pdf(fitRear.pd, xg), 'b-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Rear fit: %s', fitRear.bestName));
end

xlabel('Avg Specific Power (W/cm^2)'); ylabel('Probability Density');
title(titleStr, 'Interpreter', 'none');
legend('Location', 'best', 'FontSize', 7);
end


function row = format_stats_row(groupLabel, cornerLabel, fitResult)
% Formats one fit_distribution_summary result into a uitable row.
if isempty(fitResult.bestParams)
    paramsStr = '-';
else
    paramsStr = sprintf('%.4g, ', fitResult.bestParams);
    paramsStr = paramsStr(1:end-2);   % drop trailing ", "
end
row = {groupLabel, cornerLabel, fitResult.n, ...
    round(fitResult.mean, 3), round(fitResult.median, 3), round(fitResult.std, 3), ...
    round(fitResult.skew, 3), fitResult.bestName, paramsStr};
end


function modelResult = run_model_analysis(datasets, model, modelIdx, C)
% Runs the full per-model analysis pipeline (simulation -> events -> all
% figures) for ONE padFrac model, using the already-loaded/parsed/
% event-detected per-file data cached in `datasets` (see "LOAD & CACHE
% PER-FILE DATA" in the main script). Every figure/title/console header
% is tagged with modelIdx/model.name so multiple models' outputs stay
% distinguishable when run back to back. Returns a compact modelResult
% (name + per-event avg_q vectors) used to build the final cross-model
% summary figure after all models have been run.

tag = sprintf('[M%d: %s]', modelIdx, model.name);
fprintf('\n================ %s ================\n', tag);

VehicleMass = C.VehicleMass; RotorMass_front = C.RotorMass_front; RotorMass_rear = C.RotorMass_rear;
RotorArea_front = C.RotorArea_front; RotorArea_rear = C.RotorArea_rear; I = C.I;
A_pad_front_cm2 = C.A_pad_front_cm2; A_pad_rear_cm2 = C.A_pad_rear_cm2;
pads_per_corner = C.pads_per_corner;
BrakeFrac = C.BrakeFrac; CalibrationFactor = C.CalibrationFactor; min_pressure = C.min_pressure;
min_omega_wheel_rad_s = C.min_omega_wheel_rad_s; min_pressure_muratio_psi = C.min_pressure_muratio_psi;
mu_ratio_fade_threshold = C.mu_ratio_fade_threshold; n_temp_bins = C.n_temp_bins; n_flux_bins = C.n_flux_bins;
min_samples_per_bin = C.min_samples_per_bin; lagged_flux_window_s = C.lagged_flux_window_s;
RegenFailure_Speed_mps = C.RegenFailure_Speed_mps; RegenFailure_Decel_g = C.RegenFailure_Decel_g;
RegenFailure_TbiasFront = C.RegenFailure_TbiasFront;
maxSpeedSeen = C.maxSpeedSeen; sumTbias = C.sumTbias; countTbias = C.countTbias;

events_front = struct('dataset', {}, 't_rel', {}, 'q_inst', {}, 'avg_q', {}, 't_decel', {}, ...
    'driveType', {}, 'start_speed_mph', {}, 'decel_g', {}, 'regen_pct', {}, 'T_rotor_start_F', {});
events_rear  = struct('dataset', {}, 't_rel', {}, 'q_inst', {}, 'avg_q', {}, 't_decel', {}, ...
    'driveType', {}, 'start_speed_mph', {}, 'decel_g', {}, 'regen_pct', {}, 'T_rotor_start_F', {});

fade_samples_front = zeros(0, 4);
fade_samples_rear  = zeros(0, 4);

for k = 1:numel(datasets)
    ds = datasets(k);
    derived = ds.derived;
    t = ds.t;
    Edrag = ds.Edrag;

    % Per-corner wheel speeds [fl fr rl rr] for the strictly per-corner
    % rotational-KE sum inside simulate_pad_power.
    omega_all = [derived.fl_omega_wheel(:), derived.fr_omega_wheel(:), ...
                 derived.rl_omega_wheel(:), derived.rr_omega_wheel(:)];
    % Each dataset's own ambient / measured start temperature, matching the
    % conditions BrakeCoeffOptimizer fit these coefficients under.
    TambK_ds = ds.TambK;
    Tstart_front_K = (derived.fr_temp_F(1) - 32) * (5/9) + 273.15;
    Tstart_rear_K  = (derived.rr_temp_F(1) - 32) * (5/9) + 273.15;

    sim_front = simulate_pad_power(t, derived.velx, derived.frontpressure, derived.Tbias_brake, ...
        model.x1f, model.b1f, model.fun, derived.front_axle_regen_power, Edrag, ...
        derived.fl_omega_wheel, derived.fr_omega_wheel, omega_all, ...
        VehicleMass, RotorMass_front, RotorArea_front, I, TambK_ds, ...
        A_pad_front_cm2, pads_per_corner, BrakeFrac, CalibrationFactor, min_pressure, 0.5, Tstart_front_K, ...
        derived.T_predicted_front, min_omega_wheel_rad_s, min_pressure_muratio_psi);

    sim_rear = simulate_pad_power(t, derived.velx, derived.rearpressure, 1 - derived.Tbias_brake, ...
        model.x1r, model.b1r, model.fun, derived.rear_axle_regen_power, Edrag, ...
        derived.rl_omega_wheel, derived.rr_omega_wheel, omega_all, ...
        VehicleMass, RotorMass_rear, RotorArea_rear, I, TambK_ds, ...
        A_pad_rear_cm2, pads_per_corner, BrakeFrac, CalibrationFactor, min_pressure, 0.5, Tstart_rear_K, ...
        derived.T_predicted_rear, min_omega_wheel_rad_s, min_pressure_muratio_psi);

    dt_ds = median(diff(t), 'omitnan');
    lag_window_samples = max(1, round(lagged_flux_window_s / max(dt_ds, eps)));
    q_lag_front = movmean(sim_front.q_inst, [lag_window_samples-1, 0]);
    q_lag_rear  = movmean(sim_rear.q_inst,  [lag_window_samples-1, 0]);

    event_ranges = ds.event_ranges;

    for e = 1:size(event_ranges, 1)
        i0 = event_ranges(e, 1);
        i1 = event_ranges(e, 2);
        t_decel = t(i1) - t(i0);
        if t_decel <= 0
            continue
        end
        idx = (i0+1):i1;   % energies are defined on steps 2:end

        start_speed_mps = derived.velx(i0);
        end_speed_mps   = derived.velx(i1);
        decel_g         = (start_speed_mps - end_speed_mps) / t_decel / 9.81;
        start_speed_mph = start_speed_mps * 2.23694;

        t_rel_f = t(idx) - t(i0);
        q_inst_f = sim_front.q_inst(idx);
        total_E_f = sum(sim_front.pad_energy(idx));
        avg_q_f = total_E_f / t_decel / A_pad_front_cm2;
        regen_pct_f = 100 * sum(sim_front.regen_energy(idx)) / ...
            max(sum(sim_front.regen_energy(idx)) + sum(sim_front.friction_energy(idx)), eps);

        t_rel_r = t(idx) - t(i0);
        q_inst_r = sim_rear.q_inst(idx);
        total_E_r = sum(sim_rear.pad_energy(idx));
        avg_q_r = total_E_r / t_decel / A_pad_rear_cm2;
        regen_pct_r = 100 * sum(sim_rear.regen_energy(idx)) / ...
            max(sum(sim_rear.regen_energy(idx)) + sum(sim_rear.friction_energy(idx)), eps);

        events_front(end+1) = struct('dataset', ds.name, 't_rel', t_rel_f, ...
            'q_inst', q_inst_f, 'avg_q', avg_q_f, 't_decel', t_decel, 'driveType', ds.driveType, ...
            'start_speed_mph', start_speed_mph, 'decel_g', decel_g, 'regen_pct', regen_pct_f, ...
            'T_rotor_start_F', derived.fr_temp_F(i0)); %#ok<AGROW>
        events_rear(end+1) = struct('dataset', ds.name, 't_rel', t_rel_r, ...
            'q_inst', q_inst_r, 'avg_q', avg_q_r, 't_decel', t_decel, 'driveType', ds.driveType, ...
            'start_speed_mph', start_speed_mph, 'decel_g', decel_g, 'regen_pct', regen_pct_r, ...
            'T_rotor_start_F', derived.rr_temp_F(i0)); %#ok<AGROW>

        % Pool per-sample fade-analysis data. NaN entries (gated out by
        % the omega/pressure thresholds inside simulate_pad_power) are
        % dropped here rather than propagated into the binning step.
        mu_f = sim_front.mu_ratio(idx);
        keep_f = isfinite(mu_f);
        fade_samples_front = [fade_samples_front; ...
            q_inst_f(keep_f), q_lag_front(idx(keep_f)), derived.fr_temp_F(idx(keep_f)), mu_f(keep_f)]; %#ok<AGROW>

        mu_r = sim_rear.mu_ratio(idx);
        keep_r = isfinite(mu_r);
        fade_samples_rear = [fade_samples_rear; ...
            q_inst_r(keep_r), q_lag_rear(idx(keep_r)), derived.rr_temp_F(idx(keep_r)), mu_r(keep_r)]; %#ok<AGROW>
    end
end

nEvents = numel(events_front);
if nEvents == 0
    warning('%s: no braking events found - skipping figures for this model.', tag);
    modelResult = struct('name', model.name, 'avg_q_front', [], 'avg_q_rear', []);
    return
end

%% ---- OUTPUT 1: SPECIFIC POWER vs. TIME-SINCE-START ----
figure('Name', sprintf('%s Specific Power vs. Time Since Braking Event Start', tag));

subplot(2,1,1); hold on; grid on;
for e = 1:nEvents
    plot(events_front(e).t_rel, events_front(e).q_inst, 'r-', 'LineWidth', 0.75, 'Color', [1 0 0 0.35]);
end
xlabel('Time since event start (s)'); ylabel('Specific Power (W/cm^2)');
title(sprintf('%s Front Pad Specific Power - %d events', tag, nEvents), 'Interpreter', 'none');

subplot(2,1,2); hold on; grid on;
for e = 1:nEvents
    plot(events_rear(e).t_rel, events_rear(e).q_inst, 'b-', 'LineWidth', 0.75, 'Color', [0 0 1 0.35]);
end
xlabel('Time since event start (s)'); ylabel('Specific Power (W/cm^2)');
title(sprintf('%s Rear Pad Specific Power - %d events', tag, nEvents), 'Interpreter', 'none');

%% ---- OUTPUT 2: AVERAGE SPECIFIC POWER ----
avg_q_front_all = mean([events_front.avg_q]);
avg_q_rear_all  = mean([events_rear.avg_q]);
avg_q_combined  = mean([[events_front.avg_q], [events_rear.avg_q]]);

fprintf('\n%s SPECIFIC POWER SUMMARY\n', tag);
fprintf('Events analyzed: %d\n', nEvents);
fprintf('Average specific power, FRONT pads: %.3f W/cm^2\n', avg_q_front_all);
fprintf('Average specific power, REAR pads:  %.3f W/cm^2\n', avg_q_rear_all);
fprintf('Average specific power, COMBINED (front+rear pooled, same pad material): %.3f W/cm^2\n', avg_q_combined);
fprintf('Peak per-event average, FRONT: %.3f W/cm^2\n', max([events_front.avg_q]));
fprintf('Peak per-event average, REAR:  %.3f W/cm^2\n', max([events_rear.avg_q]));
fprintf('Peak instantaneous, FRONT: %.3f W/cm^2\n', max(cellfun(@(x) max([x; 0]), {events_front.q_inst})));
fprintf('Peak instantaneous, REAR:  %.3f W/cm^2\n', max(cellfun(@(x) max([x; 0]), {events_rear.q_inst})));

%% ---- OUTPUT 3: BOX-AND-WHISKER PLOT ----
front_vals = [events_front.avg_q]';
rear_vals  = [events_rear.avg_q]';
all_vals   = [front_vals; rear_vals];
grp        = [repmat({'Front'}, numel(front_vals), 1); repmat({'Rear'}, numel(rear_vals), 1)];

figure('Name', sprintf('%s Average Specific Power per Braking Event', tag));
boxplot(all_vals, grp);
ylabel('Average Specific Power (W/cm^2)');
title(sprintf('%s Per-Event Average Specific Power (n = %d events)', tag, nEvents), 'Interpreter', 'none');
grid on;

%% ---- OUTPUT 4: FADE CORRELATION HEATMAPS ----
% mu_ratio vs specific power, with temperature visible (not collapsed
% away) - mu has its own normal temperature dependence independent of
% fade, so a plain 2D scatter of mu_ratio vs flux alone would be
% misleading. Two flux definitions (instantaneous and lagged) are shown
% side by side per corner.
figure('Name', sprintf('%s Fade Correlation: mu Ratio vs Specific Power and Temperature', tag));

subplot(2,2,1);
plotFadeHeatmap(fade_samples_front(:,1), fade_samples_front(:,3), fade_samples_front(:,4), ...
    sprintf('%s Front - Instantaneous Flux', tag));
subplot(2,2,2);
plotFadeHeatmap(fade_samples_front(:,2), fade_samples_front(:,3), fade_samples_front(:,4), ...
    sprintf('%s Front - %.1fs Lagged Flux', tag, lagged_flux_window_s));
subplot(2,2,3);
plotFadeHeatmap(fade_samples_rear(:,1), fade_samples_rear(:,3), fade_samples_rear(:,4), ...
    sprintf('%s Rear - Instantaneous Flux', tag));
subplot(2,2,4);
plotFadeHeatmap(fade_samples_rear(:,2), fade_samples_rear(:,3), fade_samples_rear(:,4), ...
    sprintf('%s Rear - %.1fs Lagged Flux', tag, lagged_flux_window_s));

%% ---- OUTPUT 5: FADE-ONSET THRESHOLD DETECTION ----
% mu_ratio = mu_actual/mu_nominal, where mu_actual comes from dynamics
% (wheel deceleration), NOT from pressure+assumed mu - so this ratio
% dropping below 1.0 is a real friction-coefficient shortfall, not a
% modeling artifact. Binned by temperature first (to control for mu's
% normal temperature dependence), then by specific power within each
% temperature bin, to find where the ratio drops below
% mu_ratio_fade_threshold and stays there as flux keeps increasing.
fadeFront = analyzeFadeOnset(fade_samples_front, mu_ratio_fade_threshold, ...
    n_temp_bins, n_flux_bins, min_samples_per_bin);
fadeRear  = analyzeFadeOnset(fade_samples_rear, mu_ratio_fade_threshold, ...
    n_temp_bins, n_flux_bins, min_samples_per_bin);

fprintf('\n%s FADE-ONSET ANALYSIS\n', tag);
fprintf('Fade defined as mu_actual/mu_nominal dropping below %.2f and staying there as flux increases.\n', ...
    mu_ratio_fade_threshold);
fprintf('Computed within %d temperature bins to control for mu''s own normal temperature dependence.\n\n', ...
    n_temp_bins);

report_fade_corner('FRONT', fadeFront, A_pad_front_cm2, avg_q_front_all, ...
    max([events_front.avg_q]), lagged_flux_window_s);
report_fade_corner('REAR', fadeRear, A_pad_rear_cm2, avg_q_rear_all, ...
    max([events_rear.avg_q]), lagged_flux_window_s);

%% ---- OUTPUT 6: DRIVE-TYPE BREAKDOWN - PROBABILITY DISTRIBUTIONS + STATS TABLE ----
driveCats = {'Endurance', 'Autocross', 'Other/Misc'};
frontDriveTypes = {events_front.driveType};
rearDriveTypes  = {events_rear.driveType};
frontAvgQ = [events_front.avg_q];
rearAvgQ  = [events_rear.avg_q];

figure('Name', sprintf('%s Avg Specific Power by Drive Type - Distributions', tag));
tableRows = {};
for c = 1:numel(driveCats)
    maskF = strcmp(frontDriveTypes, driveCats{c});
    maskR = strcmp(rearDriveTypes, driveCats{c});

    subplot(2,3,c);
    [fitF, fitR] = plot_distribution_overlay(frontAvgQ(maskF), rearAvgQ(maskR), ...
        sprintf('%s\n(n_{front}=%d, n_{rear}=%d)', driveCats{c}, sum(maskF), sum(maskR)));

    tableRows(end+1,:) = format_stats_row(driveCats{c}, 'Front', fitF); %#ok<AGROW>
    tableRows(end+1,:) = format_stats_row(driveCats{c}, 'Rear',  fitR); %#ok<AGROW>
end
% Bottom row of the 2x3 grid is deliberately left empty (only subplot
% positions 1-3 are ever used above) so the table can occupy that space
% without fighting the axes for room.
uitable('Parent', gcf, 'Units', 'normalized', 'Position', [0.03 0.03 0.94 0.40], ...
    'Data', tableRows, ...
    'ColumnName', {'Drive Type', 'Corner', 'n', 'Mean', 'Median', 'Std', 'Skew', 'Best Fit', 'Params'}, ...
    'ColumnWidth', {80, 50, 40, 70, 70, 70, 60, 90, 220}, 'RowName', []);

fprintf('\n%s DRIVE-TYPE BREAKDOWN\n', tag);
for c = 1:numel(driveCats)
    maskF = strcmp(frontDriveTypes, driveCats{c});
    maskR = strcmp(rearDriveTypes, driveCats{c});
    fprintf('%-12s Front: n=%d   Rear: n=%d\n', driveCats{c}, sum(maskF), sum(maskR));
end
fprintf('\n');

%% ---- OUTPUT 7: SPECIFIC POWER CORRELATIONS ----
figure('Name', sprintf('%s Specific Power Correlations', tag));
corrSpecs = {
    'decel_g',         'Deceleration (g)'
    'start_speed_mph', 'Starting Speed (mph)'
    'regen_pct',       'Regen % of Braking Energy'
    't_decel',         'Event Duration (s)'
    'T_rotor_start_F', 'Rotor Temp at Event Start (deg F)'
    };
for p = 1:size(corrSpecs, 1)
    fieldName = corrSpecs{p,1};
    subplot(2,3,p); hold on; grid on;
    scatter([events_front.(fieldName)], frontAvgQ, 18, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.5, 'DisplayName', 'Front');
    scatter([events_rear.(fieldName)], rearAvgQ, 18, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.5, 'DisplayName', 'Rear');
    xlabel(corrSpecs{p,2}); ylabel('Specific Power (W/cm^2)');
    title(sprintf('%s vs. Specific Power', corrSpecs{p,2}), 'Interpreter', 'none');
    if p == 1
        legend('Location', 'best');
    end
end
sgtitle(tag, 'Interpreter', 'none');

%% ---- OUTPUT 8: REGEN-FAILURE WORST-CASE MARGIN ----
fprintf('\n%s REGEN-FAILURE WORST-CASE MARGIN\n', tag);
fprintf(['NOTE: the assumptions below are editable placeholders (see ', ...
    '"REGEN-FAILURE WORST-CASE ASSUMPTIONS" near the top of the script) - ', ...
    'review before trusting the back-calculated pad area.\n\n']);

if isempty(RegenFailure_Speed_mps)
    v0 = maxSpeedSeen;
    fprintf('Speed assumption: max speed observed in loaded data = %.1f m/s (%.1f mph)\n', v0, v0*2.23694);
else
    v0 = RegenFailure_Speed_mps;
    fprintf('Speed assumption: user-specified = %.1f m/s (%.1f mph)\n', v0, v0*2.23694);
end

if isempty(RegenFailure_TbiasFront)
    if countTbias > 0
        TbiasF = sumTbias / countTbias;
    else
        TbiasF = 0.5;
        warning('No valid Tbias_brake samples observed - defaulting regen-failure split to 0.5/0.5.');
    end
    fprintf('Front/rear split: mean observed Tbias_brake = %.3f\n', TbiasF);
else
    TbiasF = RegenFailure_TbiasFront;
    fprintf('Front/rear split: user-specified Tbias = %.3f\n', TbiasF);
end

a_decel = RegenFailure_Decel_g * 9.81;
fprintf('Assumed constant deceleration: %.2f g (%.2f m/s^2)\n', RegenFailure_Decel_g, a_decel);

F_total = VehicleMass * a_decel;    % N, ALL via friction (zero regen assumed)
P_total_peak = F_total * v0;        % W, peak power at t=0 (highest speed, constant-decel stop)
% Divided all the way down to ONE PAD so it can be compared against qOnset,
% which is a per-single-pad flux. (x TbiasF -> one axle, x0.5 -> one corner,
% / pads_per_corner -> one pad.) These two MUST stay on the same basis: the
% margin factor below is P_onset/P_peak with pad area cancelling, so a
% mismatch here silently scales the reported margin.
P_front_peak = P_total_peak * TbiasF       * 0.5 / pads_per_corner;   % per SINGLE front pad
P_rear_peak  = P_total_peak * (1 - TbiasF) * 0.5 / pads_per_corner;   % per SINGLE rear pad

fprintf('Peak total friction power at t=0: %.1f kW\n', P_total_peak/1000);
fprintf('Peak per-PAD power: front = %.2f kW, rear = %.2f kW\n\n', P_front_peak/1000, P_rear_peak/1000);

qOnsetFront = fadeFront.q_inst.onsetOverall;
qOnsetRear  = fadeRear.q_inst.onsetOverall;

if isfinite(qOnsetFront) && qOnsetFront > 0
    A_required_front_cm2 = P_front_peak / qOnsetFront;
    fprintf(['Front: required pad area to stay at/under fade-onset flux = %.1f cm^2 ', ...
        '(currently %.1f cm^2, margin factor %.2fx)\n'], ...
        A_required_front_cm2, A_pad_front_cm2, A_pad_front_cm2/A_required_front_cm2);
else
    fprintf('Front: fade-onset flux not determined from the loaded data - cannot back-calculate required area.\n');
end

if isfinite(qOnsetRear) && qOnsetRear > 0
    A_required_rear_cm2 = P_rear_peak / qOnsetRear;
    fprintf(['Rear: required pad area to stay at/under fade-onset flux = %.1f cm^2 ', ...
        '(currently %.1f cm^2, margin factor %.2fx)\n'], ...
        A_required_rear_cm2, A_pad_rear_cm2, A_pad_rear_cm2/A_required_rear_cm2);
else
    fprintf('Rear: fade-onset flux not determined from the loaded data - cannot back-calculate required area.\n');
end

modelResult = struct('name', model.name, 'avg_q_front', frontAvgQ, 'avg_q_rear', rearAvgQ);
end
