%% ================================================================
%  Brake Pad Specific Power (Heat Flux) & Fade-Onset Analysis
% ================================================================
%  Uses FIXED, already-calibrated rotor cooling and pad-fraction
%  coefficients (set as constants below - see "CALIBRATED MODEL
%  CONSTANTS") to walk one or more driving datasets and compute, for
%  every identified braking event:
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
%  Outputs:
%    1) Specific power (W/cm^2) vs. time-since-event-start, for every
%       braking event found, front and rear pads plotted separately.
%    2) The average specific power across all events (front & rear).
%    3) A box-and-whisker plot of each event's average specific power
%       (front vs. rear).
%    4) Heatmaps of mu_actual/mu_nominal vs. specific power and
%       temperature, for both instantaneous and lagged flux.
%    5) The empirically-derived fade-onset specific power, front/rear.
%    6) Average specific power broken out by drive type (endurance,
%       autocross, other/misc), auto-detected from preprocessed file
%       headers where available.
%    7) A regen-braking-failure worst-case margin check: back-calculates
%       the pad area needed to stay under the fade-onset flux.
%
%  Datasets are selected via a multi-select file picker, same
%  21/22/27-column auto-detected formats as BrakeDataAnalysis.m.
% ================================================================

clc; clear; close all

%% ================== CALIBRATED MODEL CONSTANTS ==================
% >>> Replace these with your finalized fit from BrakeCoeffOptimizer.m <<<

% Rotor cooling coefficient: h_w(v) = x1*v + b1   [W/m^2-K]
x1 = 2.0006;
b1 = 40.7137;

% Pad energy fraction: PadFrac = f(T_rotor [K], P_applied [psi])
% Default below matches the "linear in T only" baseline model.
% If you settled on a pressure-dependent model instead, edit padfrac_fun
% and the coefficients to match (see commented examples).
x2 = 0.001005;    % 1/K
b2 = -0.5385;      % intercept
padfrac_fun = @(T, P) x2.*T + b2;

% --- Example alternates (uncomment + set coefficients as needed) ---
% x2 = ...; x3 = ...; b2 = ...;
% padfrac_fun = @(T, P) x2.*T + x3.*P + b2;                  % linear, indep. T & P
%
% x2 = ...; x3 = ...; x4 = ...; b2 = ...;
% padfrac_fun = @(T, P) x2.*T + x3.*P + x4.*T.*P + b2;       % T*P interaction
%
% x2 = ...; x2q = ...; x3 = ...; b2 = ...;
% padfrac_fun = @(T, P) x2.*T + x2q.*T.^2 + x3.*P + b2;      % quadratic in T

%% ================== VEHICLE / BRAKE CONSTANTS ==================
VehicleMass     = 259;      % kg
RotorMass_front = 0.5107;   % kg
RotorMass_rear  = 0.3343;   % kg
RotorArea_front = 0.038;    % m^2 (rotor face area, for cooling calc)
RotorArea_rear  = 0.0226;   % m^2
I               = 0.30754;  % rotational inertia kg*m^2
WheelR          = 0.213;    % meters
gear_ratio      = 12.97;
TambC           = 25;       % degC
BrakeFrac       = 0.82;     % mechanical brake fraction (matches brake_temp_sim.m)
CalibrationFactor = 0.8;    % matches the fixed 0.8 factor in brake_temp_sim.m

% Pad areas -- PER SINGLE PAD (not combined per axle), from spec sheet
A_pad_front_mm2 = 1520;
A_pad_rear_mm2  = 650;
A_pad_front_cm2 = A_pad_front_mm2 / 100;
A_pad_rear_cm2  = A_pad_rear_mm2  / 100;

% Brake geometry (needed to reconstruct Tbias_brake from raw pressures)
front_piston_count = 6;  front_piston_dia = 0.0157226;  front_rotor_dia = 0.18542;
rear_piston_count  = 4;  rear_piston_dia  = 0.0141732;  rear_rotor_dia  = 0.18796;

mu_temp_table = [100, 200, 300, 400, 500, 600, 700, 800, 900, 950, 1100, 1200];  % degF
mu_table      = [0.45, 0.46, 0.49, 0.53, 0.55, 0.56, 0.565, 0.565, 0.57, 0.57, 0.535, 0.535];

aero_open_a   =  0.0793024;   aero_open_b   =  1.46483;    aero_open_c   = -37.23993;
aero_closed_a =  0.136247;    aero_closed_b =  2.53449;    aero_closed_c = -67.44975;

velx_threshold = 34;  % m/s, cleaning threshold
min_pressure   = 5;   % psi, threshold to consider a corner's brake engaged

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

aeroParams = struct('front_piston_count', front_piston_count, 'front_piston_dia', front_piston_dia, ...
    'rear_piston_count', rear_piston_count, 'rear_piston_dia', rear_piston_dia, ...
    'front_rotor_dia', front_rotor_dia, 'rear_rotor_dia', rear_rotor_dia, ...
    'mu_temp_table', mu_temp_table, 'mu_table', mu_table, ...
    'aero_open_a', aero_open_a, 'aero_open_b', aero_open_b, 'aero_open_c', aero_open_c, ...
    'aero_closed_a', aero_closed_a, 'aero_closed_b', aero_closed_b, 'aero_closed_c', aero_closed_c);
gearParams = struct('gear_ratio', gear_ratio, 'velx_threshold', velx_threshold);

TambK = TambC + 273.15;

% Pooled event records across all datasets
events_front = struct('dataset', {}, 't_rel', {}, 'q_inst', {}, 'avg_q', {}, 't_decel', {}, 'driveType', {});
events_rear  = struct('dataset', {}, 't_rel', {}, 'q_inst', {}, 'avg_q', {}, 't_decel', {}, 'driveType', {});

% Pooled per-SAMPLE fade-analysis data (finer than per-event averages -
% every active-braking, omega-gated sample across every event/dataset).
% Columns: [q_inst, q_lag, T_rotor_F, mu_ratio]
fade_samples_front = zeros(0, 4);
fade_samples_rear  = zeros(0, 4);

% Running accumulators used as defaults for the regen-failure assumptions
% (max speed observed, and mean Tbias_brake across all loaded data).
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

    t = parsed.t;

    %% ---- Run the full-timeseries simulation (front & rear) ----
    sim_front = simulate_pad_power(t, derived.velx, derived.frontpressure, derived.Tbias_brake, ...
        x1, b1, padfrac_fun, derived.total_regen_power, Edrag, ...
        derived.fl_omega_wheel, derived.fr_omega_wheel, ...
        VehicleMass, RotorMass_front, RotorArea_front, I, WheelR, TambK, ...
        A_pad_front_cm2, BrakeFrac, CalibrationFactor, min_pressure, 0.5, TambC, ...
        derived.T_predicted_front, min_omega_wheel_rad_s);

    sim_rear = simulate_pad_power(t, derived.velx, derived.rearpressure, 1 - derived.Tbias_brake, ...
        x1, b1, padfrac_fun, derived.total_regen_power, Edrag, ...
        derived.rl_omega_wheel, derived.rr_omega_wheel, ...
        VehicleMass, RotorMass_rear, RotorArea_rear, I, WheelR, TambK, ...
        A_pad_rear_cm2, BrakeFrac, CalibrationFactor, min_pressure, 0.5, TambC, ...
        derived.T_predicted_rear, min_omega_wheel_rad_s);

    % Lagged (trailing moving-average) specific power - candidate flux
    % definition for the fade-onset correlation, compared against the
    % instantaneous q_inst since surface-layer formation isn't
    % necessarily instantaneous.
    dt_ds = median(diff(t), 'omitnan');
    lag_window_samples = max(1, round(lagged_flux_window_s / max(dt_ds, eps)));
    q_lag_front = movmean(sim_front.q_inst, [lag_window_samples-1, 0]);
    q_lag_rear  = movmean(sim_rear.q_inst,  [lag_window_samples-1, 0]);

    %% ---- Determine this dataset's drive type (endurance/autocross/other) ----
    driveType = detectDriveType(fname);
    fprintf('  Drive type: %s\n', driveType);

    maxSpeedSeen = max(maxSpeedSeen, max(derived.velx, [], 'omitnan'));
    validTbias = derived.Tbias_brake(isfinite(derived.Tbias_brake) & derived.Tbias_brake > 0);
    sumTbias = sumTbias + sum(validTbias);
    countTbias = countTbias + numel(validTbias);

    %% ---- Identify vehicle-level braking events ----
    DS = [0; diff(derived.velx)];
    DS(DS < -2.5) = 0;  % matches the >2.5g exclusion used in the sim
    active = (derived.frontpressure > min_pressure | derived.rearpressure > min_pressure) & (DS < 0);

    event_ranges = find_events(active, t, bridge_gap_s, min_event_dur_s);
    fprintf('  Found %d braking event(s).\n', size(event_ranges, 1));

    for e = 1:size(event_ranges, 1)
        i0 = event_ranges(e, 1);
        i1 = event_ranges(e, 2);
        t_decel = t(i1) - t(i0);
        if t_decel <= 0
            continue
        end
        idx = (i0+1):i1;   % energies are defined on steps 2:end

        t_rel_f = t(idx) - t(i0);
        q_inst_f = sim_front.q_inst(idx);
        total_E_f = sum(sim_front.pad_energy(idx));
        avg_q_f = total_E_f / t_decel / A_pad_front_cm2;

        t_rel_r = t(idx) - t(i0);
        q_inst_r = sim_rear.q_inst(idx);
        total_E_r = sum(sim_rear.pad_energy(idx));
        avg_q_r = total_E_r / t_decel / A_pad_rear_cm2;

        events_front(end+1) = struct('dataset', files{k}, 't_rel', t_rel_f, ...
            'q_inst', q_inst_f, 'avg_q', avg_q_f, 't_decel', t_decel, 'driveType', driveType); %#ok<SAGROW>
        events_rear(end+1) = struct('dataset', files{k}, 't_rel', t_rel_r, ...
            'q_inst', q_inst_r, 'avg_q', avg_q_r, 't_decel', t_decel, 'driveType', driveType); %#ok<SAGROW>

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
    error('No braking events found across the selected dataset(s).');
end

%% ================== OUTPUT 1: SPECIFIC POWER vs. TIME-SINCE-START ==================
figure('Name', 'Specific Power vs. Time Since Braking Event Start');

subplot(2,1,1); hold on; grid on;
for e = 1:nEvents
    plot(events_front(e).t_rel, events_front(e).q_inst, 'r-', 'LineWidth', 0.75, 'Color', [1 0 0 0.35]);
end
xlabel('Time since event start (s)'); ylabel('Specific Power (W/cm^2)');
title(sprintf('Front Pad Specific Power - %d events', nEvents));

subplot(2,1,2); hold on; grid on;
for e = 1:nEvents
    plot(events_rear(e).t_rel, events_rear(e).q_inst, 'b-', 'LineWidth', 0.75, 'Color', [0 0 1 0.35]);
end
xlabel('Time since event start (s)'); ylabel('Specific Power (W/cm^2)');
title(sprintf('Rear Pad Specific Power - %d events', nEvents));

%% ================== OUTPUT 2: AVERAGE SPECIFIC POWER ==================
avg_q_front_all = mean([events_front.avg_q]);
avg_q_rear_all  = mean([events_rear.avg_q]);
avg_q_combined  = mean([[events_front.avg_q], [events_rear.avg_q]]);

fprintf('\n================ SPECIFIC POWER SUMMARY ================\n');
fprintf('Events analyzed: %d\n', nEvents);
fprintf('Average specific power, FRONT pads: %.3f W/cm^2\n', avg_q_front_all);
fprintf('Average specific power, REAR pads:  %.3f W/cm^2\n', avg_q_rear_all);
fprintf('Average specific power, COMBINED (front+rear pooled, same pad material): %.3f W/cm^2\n', avg_q_combined);
fprintf('Peak per-event average, FRONT: %.3f W/cm^2\n', max([events_front.avg_q]));
fprintf('Peak per-event average, REAR:  %.3f W/cm^2\n', max([events_rear.avg_q]));
fprintf('Peak instantaneous, FRONT: %.3f W/cm^2\n', max(cellfun(@(x) max([x; 0]), {events_front.q_inst})));
fprintf('Peak instantaneous, REAR:  %.3f W/cm^2\n', max(cellfun(@(x) max([x; 0]), {events_rear.q_inst})));

%% ================== OUTPUT 3: BOX-AND-WHISKER PLOT ==================
front_vals = [events_front.avg_q]';
rear_vals  = [events_rear.avg_q]';
all_vals   = [front_vals; rear_vals];
grp        = [repmat({'Front'}, numel(front_vals), 1); repmat({'Rear'}, numel(rear_vals), 1)];

figure('Name', 'Average Specific Power per Braking Event');
boxplot(all_vals, grp);
ylabel('Average Specific Power (W/cm^2)');
title(sprintf('Per-Event Average Specific Power (n = %d events)', nEvents));
grid on;

%% ================== OUTPUT 4: FADE CORRELATION HEATMAPS ==================
% mu_ratio vs specific power, with temperature visible (not collapsed
% away) - mu has its own normal temperature dependence independent of
% fade, so a plain 2D scatter of mu_ratio vs flux alone would be
% misleading. Two flux definitions (instantaneous and lagged) are shown
% side by side per corner.
figure('Name', 'Fade Correlation: mu Ratio vs Specific Power and Temperature');

subplot(2,2,1);
plotFadeHeatmap(fade_samples_front(:,1), fade_samples_front(:,3), fade_samples_front(:,4), ...
    'Front - Instantaneous Flux');
subplot(2,2,2);
plotFadeHeatmap(fade_samples_front(:,2), fade_samples_front(:,3), fade_samples_front(:,4), ...
    sprintf('Front - %.1fs Lagged Flux', lagged_flux_window_s));
subplot(2,2,3);
plotFadeHeatmap(fade_samples_rear(:,1), fade_samples_rear(:,3), fade_samples_rear(:,4), ...
    'Rear - Instantaneous Flux');
subplot(2,2,4);
plotFadeHeatmap(fade_samples_rear(:,2), fade_samples_rear(:,3), fade_samples_rear(:,4), ...
    sprintf('Rear - %.1fs Lagged Flux', lagged_flux_window_s));

%% ================== OUTPUT 5: FADE-ONSET THRESHOLD DETECTION ==================
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

fprintf('\n================ FADE-ONSET ANALYSIS ================\n');
fprintf('Fade defined as mu_actual/mu_nominal dropping below %.2f and staying there as flux increases.\n', ...
    mu_ratio_fade_threshold);
fprintf('Computed within %d temperature bins to control for mu''s own normal temperature dependence.\n\n', ...
    n_temp_bins);

report_fade_corner('FRONT', fadeFront, A_pad_front_cm2, avg_q_front_all, ...
    max([events_front.avg_q]), lagged_flux_window_s);
report_fade_corner('REAR', fadeRear, A_pad_rear_cm2, avg_q_rear_all, ...
    max([events_rear.avg_q]), lagged_flux_window_s);

%% ================== OUTPUT 6: DRIVE-TYPE BREAKDOWN ==================
driveCats = {'Endurance', 'Autocross', 'Other/Misc'};
meanQ_front = nan(numel(driveCats),1); semQ_front = nan(numel(driveCats),1); nQ_front = zeros(numel(driveCats),1);
meanQ_rear  = nan(numel(driveCats),1); semQ_rear  = nan(numel(driveCats),1); nQ_rear  = zeros(numel(driveCats),1);

frontDriveTypes = {events_front.driveType};
rearDriveTypes  = {events_rear.driveType};
frontAvgQ = [events_front.avg_q];
rearAvgQ  = [events_rear.avg_q];

for c = 1:numel(driveCats)
    maskF = strcmp(frontDriveTypes, driveCats{c});
    maskR = strcmp(rearDriveTypes, driveCats{c});
    nQ_front(c) = sum(maskF);
    nQ_rear(c)  = sum(maskR);
    if nQ_front(c) > 0
        meanQ_front(c) = mean(frontAvgQ(maskF));
        semQ_front(c)  = std(frontAvgQ(maskF)) / sqrt(nQ_front(c));
    end
    if nQ_rear(c) > 0
        meanQ_rear(c) = mean(rearAvgQ(maskR));
        semQ_rear(c)  = std(rearAvgQ(maskR)) / sqrt(nQ_rear(c));
    end
end

figure('Name', 'Average Specific Power by Drive Type');
barVals = [meanQ_front, meanQ_rear];
b = bar(barVals);
hold on;
errorbar(b(1).XEndPoints, meanQ_front, semQ_front, 'k.', 'LineWidth', 1);
errorbar(b(2).XEndPoints, meanQ_rear,  semQ_rear,  'k.', 'LineWidth', 1);
set(gca, 'XTick', 1:numel(driveCats), 'XTickLabel', driveCats);
ylabel('Average Specific Power (W/cm^2)');
legend({'Front', 'Rear'}, 'Location', 'best');
title('Average Specific Power by Drive Type (error bars = SEM)');
grid on;

fprintf('================ DRIVE-TYPE BREAKDOWN ================\n');
for c = 1:numel(driveCats)
    fprintf('%-12s Front: %.3f W/cm^2 (n=%d events)   Rear: %.3f W/cm^2 (n=%d events)\n', ...
        driveCats{c}, meanQ_front(c), nQ_front(c), meanQ_rear(c), nQ_rear(c));
end
fprintf('\n');

%% ================== OUTPUT 7: REGEN-FAILURE WORST-CASE MARGIN ==================
fprintf('================ REGEN-FAILURE WORST-CASE MARGIN ================\n');
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
P_front_peak = P_total_peak * TbiasF * 0.5;         % per SINGLE front corner
P_rear_peak  = P_total_peak * (1 - TbiasF) * 0.5;   % per SINGLE rear corner

fprintf('Peak total friction power at t=0: %.1f kW\n', P_total_peak/1000);
fprintf('Peak per-corner power: front = %.2f kW, rear = %.2f kW\n\n', P_front_peak/1000, P_rear_peak/1000);

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

F_aero_open   = aeroParams.aero_open_a   * speed_mph.^2 + aeroParams.aero_open_b   * speed_mph + aeroParams.aero_open_c;
F_aero_closed = aeroParams.aero_closed_a * speed_mph.^2 + aeroParams.aero_closed_b * speed_mph + aeroParams.aero_closed_c;
F_aero = F_aero_open .* double(parsed.drs_state == 1) + F_aero_closed .* double(parsed.drs_state == 0);
F_aero = max(F_aero, 0);

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
% Per-corner PREDICTED (pressure x piston area x nominal mu x radius)
% torque, already computed above for Tbias_brake - exposed here for the
% dynamics-vs-predicted fade comparison. Sign is dropped (magnitude only)
% since we only need it for a ratio against the (also-positive) actual
% torque. Front/rear each have a single hydraulic circuit, so fl==fr and
% rl==rr by construction (same value both corners of an axle).
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
    total_regen_power, Edrag, omega_wheel_L, omega_wheel_R, ...
    VehicleMass, RotorMass, RotorArea, I, WheelR, TambK, A_pad_cm2, BrakeFrac, CalibrationFactor, ...
    min_pressure, corner_split, TambC, T_predicted, min_omega_wheel_rad_s) %#ok<INUSD>
% Walks the full timeseries computing rotor temp (needed for PadFrac's
% temperature dependence), per-step pad energy, and instantaneous
% specific power (W/cm^2). PadFrac uses the SAME (1-PadFrac)/PadFrac
% split convention as brake_temp_sim.m: PadFrac is the pad's own share,
% (1-PadFrac) is the rotor's share.
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

n = length(t);
RotorTempArrayK = zeros(n, 1);
RotorTempArrayK(1) = TambC + 273.15;   % assume ambient at dataset start
pad_energy = zeros(n, 1);   % J, per single pad, per step
q_inst     = zeros(n, 1);   % W/cm^2, per single pad, per step
T_actual   = nan(n, 1);     % Nm, per single corner, dynamics-based
mu_ratio   = nan(n, 1);     % T_actual / T_predicted

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
        omegaP  = (omega_wheel_L(i-1) + omega_wheel_R(i-1)) / 2;
        omegaN  = (omega_wheel_L(i)   + omega_wheel_R(i))   / 2;
        Energy2 = 4 * (0.5 * I * (omegaP^2 - omegaN^2));
        Energy  = Energy1 + Energy2;

        regen_energy    = abs(total_regen_power(i)) * tbrake;
        friction_energy = max(Energy - regen_energy, 0);
        AeroFrac  = min(Edrag(i) / max(Energy, 1), 1);

        % Rotor share (drives rotor temp forward, same as brake_temp_sim.m)
        CorrectedEnergyRotor = friction_energy * corner_split * Tbias(i) * (1 - AeroFrac) * (1 - PadFrac) * BrakeFrac * CalibrationFactor;
        deltaTK = CorrectedEnergyRotor / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = deltaTK + prevTemp;

        qout = (RotorTempArrayK(i) - TambK) / Rrotor;
        Eout = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = RotorTempArrayK(i) - deltaTKout;

        % Pad share (this is what we report as specific power)
        CorrectedEnergyPad = friction_energy * corner_split * Tbias(i) * (1 - AeroFrac) * PadFrac * BrakeFrac * CalibrationFactor;
        pad_energy(i) = CorrectedEnergyPad;
        q_inst(i)     = CorrectedEnergyPad / tbrake / A_pad_cm2;

        % Dynamics-based (mu-independent) actual torque: mechanical brake
        % power at this corner, divided by wheel angular velocity. Gated
        % on a minimum omega to avoid Power/omega blowing up near a stop.
        if omegaN > min_omega_wheel_rad_s
            T_actual_axle_power = (friction_energy / tbrake) * Tbias(i) * (1 - AeroFrac);
            T_actual(i) = T_actual_axle_power * corner_split / omegaN;
            if isfinite(T_predicted(i)) && T_predicted(i) > 1e-6
                mu_ratio(i) = T_actual(i) / T_predicted(i);
            end
        end
    else
        qout = (prevTemp - TambK) / Rrotor;
        Eout = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = prevTemp - deltaTKout;
        % pad_energy(i), q_inst(i), T_actual(i), mu_ratio(i) remain 0/NaN
    end
end

sim.RotorTempArrayK = RotorTempArrayK;
sim.pad_energy      = pad_energy;
sim.q_inst           = q_inst;
sim.T_actual         = T_actual;
sim.mu_ratio         = mu_ratio;
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
% Reads the first several lines of a data file looking for the
% "# Segment type: <kind>" header F34BrakeDataPreprocessor.m writes, and
% maps it to a coarse drive-type category. Falls back to a prompt for
% files without a recognizable header (e.g. raw/non-preprocessed input).
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
% Same griddata + contourf + imgaussfilt heatmap style used for the
% throttle/steer/yaw heatmaps, applied to mu_ratio vs (flux, temperature)
% so temperature context is never collapsed away.
valid = isfinite(flux) & isfinite(tempF) & isfinite(muRatio);
flux = flux(valid); tempF = tempF(valid); muRatio = muRatio(valid);

if numel(flux) < 10 || range(flux) <= 0 || range(tempF) <= 0
    title([titleStr ' (insufficient data)']);
    axis off;
    return
end

flux_grid = linspace(min(flux), max(flux), 50);
temp_grid = linspace(min(tempF), max(tempF), 50);
[F_grid, T_grid] = meshgrid(flux_grid, temp_grid);

Z_grid = griddata(flux, tempF, muRatio, F_grid, T_grid, 'linear');
Z_grid = imgaussfilt(Z_grid, 1.5, 'FillValues', NaN);

contourf(F_grid, T_grid, Z_grid, 20, 'LineColor', 'none');
hold on;
contour(F_grid, T_grid, Z_grid, 10, 'LineColor', 'k', 'LineWidth', 0.5);
scatter(flux, tempF, 6, 'w.', 'MarkerEdgeAlpha', 0.25);
colormap('jet');
c = colorbar; c.Label.String = 'mu_{actual} / mu_{nominal}';
xlabel('Specific Power (W/cm^2)');
ylabel('Rotor Temp (deg F, pad-face proxy)');
title(titleStr, 'Interpreter', 'none');
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