%% ================================================================
%  Regenerative Braking Energy Contribution Analysis
% ================================================================
%  Answers one question: of the energy dissipated in a braking event, how
%  much of it does regen recover?
%
%  Combines the regen physics of BrakeDataAnalysis.m (regen / mechanical /
%  aero energy accounting, torque bias, and the current-limited theoretical
%  regen torque model) with the multi-dataset front end of
%  padSpecificPower.m and BrakeCoeffOptimizer.m: a multi-select file picker,
%  21/22/27-column format auto-detection, per-file parse/derive, drive-type
%  detection (Autocross / Endurance / Other-Misc), and braking-event
%  segmentation.
%
%  A "braking event" here is ANY deceleration event - the friction brakes do
%  not have to be applied. Pure-regen stops are real braking events and
%  excluding them would bias the headline regen share downward.
%
%  Datasets recorded without regen are DROPPED entirely - they contribute to
%  no calculation and appear in no figure - so they cannot dilute the
%  averages. Sessions run without regen do exist in carData: BOTH 4-17
%  endurance sessions recover ~0.2% of their dissipated energy, against
%  10-48% for every later session. The test uses pack current where the file
%  provides it and regen energy share where it does not - see
%  "REGEN-PRESENCE SETTINGS" below for why that split is necessary.
%
%  Console output (per dataset, per drive type, then pooled):
%    - Regen contribution headline: energy-weighted AND mean-per-event
%    - Braking energy summary, aero drag during regen, regen model vs
%      actual, braking energy breakdown, DRS drag energy analysis
%
%  Figures:
%    1) Distribution of mechanical brake torque bias (all in-event samples)
%    2) Distribution of total decel torque bias (all in-event samples)
%    3) Regen torque bias vs speed
%    4) Total decel torque bias vs speed
%    5) Regen total wheel torque vs speed - one figure colored by BMS
%       current (files that have it), one colored by regen power (files
%       that do not)
%    6) Theoretical max regen torque vs speed at several current limits
%    7) Cumulative braking energy breakdown, one pie per drive type
%    8) Per-event regen / friction / aero energy distributions per drive type
%    9) Regen energy vs deceleration rate, and vs event duration
%
%  Rotor temperature is still USED (mu vs temperature drives the friction
%  brake torque) but no temperature output is produced - that lives in
%  BrakeDataAnalysis.m and padSpecificPower.m.
% ================================================================

clc; clear; close all

%% ================== VEHICLE / REGEN CONSTANTS ==================
% Everything except the current-limit sweep comes from the CENTRALIZED
% spreadsheet VehicleParameters.xlsx (see loadVehicleParams.m) so this
% script cannot drift from BrakeCoeffOptimizer.m / padSpecificPower.m.
VP = loadVehicleParams();

WheelR               = VP.WheelR;
gear_ratio           = VP.gear_ratio;
min_regen_speed      = VP.min_regen_speed;        % mph, system regen cut-off
min_energy_threshold = VP.min_energy_threshold;   % kJ, floor for % difference
current_limit_actual = VP.current_limit_actual;   % A, setting the data was recorded at
regen_power_limit    = VP.regen_power_limit;      % W, controller power limit
regen_max_curve_a    = VP.regen_max_curve_a;      % max regen torque curve T = a*mph^b
regen_max_curve_b    = VP.regen_max_curve_b;
regen_efficiency     = VP.regen_efficiency;       % utilization factor vs theoretical max

% Wheel-side total torque limit, DERIVED from the %Mn command limit - kept
% out of the spreadsheet because it is a function of gear_ratio, which is
% already stored there.
torque_limit_total_wheel = VP.torque_limit_total_MN / 100 * 9.8 * gear_ratio;  % Nm, all 4 wheels

% Current limits to compare in the theoretical-torque figure. This is a
% sweep setting rather than a vehicle property, and loadVehicleParams is
% scalar-only, so it stays here.
current_limits = [30, 35, 40];   % A
limit_colors   = {'r-', 'g-', 'b-'};

%% ================== REGEN-PRESENCE SETTINGS ==================
% Which samples count as "regenerating", and which whole datasets get
% dropped for having no regen at all, depends on whether the file carries
% pack current. It is NOT a function of column count:
%
%   - RAW logger files (carData\!F34raw) have 2_BMS_InstCurrentFilt[A].
%   - PREPROCESSED exports from F34BrakeDataPreprocessor.m DROP it - the
%     column in that position is Steering_pct. Front/rear pressure, the
%     motor channels, velocity and accel all still line up, which is why
%     the positional parser below is shared by both.
%
% So the BMS column is located BY NAME from the channel-name header row,
% and both the per-sample regen index and the dataset gate have a
% torque/energy fallback for files without it.
min_regen_current_A    = -5;   % A, pack current at or below this = genuinely charging.
                               % Not simply "< 0": the no-regen raw session
                               % (First Endurance No Aero No Regen 4-17) still
                               % dips to -3.6 A on sensor noise, while sessions
                               % with regen reach -36 to -86 A.
regen_torque_threshold = -10;  % Nm wheel-side, sample counts as regenerating below
                               % this when no BMS channel exists (the value
                               % BrakeDataAnalysis.m defines but never uses)
min_regen_share_pct    = 2;    % %, non-BMS files: regen must be at least this share
                               % of the file's dissipated energy to be kept.
                               % Sessions without regen measure 0.2-0.3%;
                               % sessions with it measure 10-48%.

%% ================== EVENT-DETECTION SETTINGS ==================
% A braking event = decelerating harder than decel_accel_threshold while
% above min_event_speed_mph. Short gaps are bridged so one stop with a
% modulated pedal isn't split into fragments, and very short blips are
% dropped - the same treatment padSpecificPower.m applies to its
% (pressure-gated) events.
decel_accel_threshold = -1.5;  % m/s^2, longitudinal accel below this = braking
min_event_speed_mph   = 5;     % mph, ignore crawl-speed noise
bridge_gap_s          = 0.2;   % bridge gaps in braking shorter than this
min_event_dur_s       = 0.2;   % drop merged events shorter than this

%% ================== FIGURE GATING SETTINGS ==================
% Display-only filters. None of these alter a computed energy or bias -
% they only decide which samples get drawn.
min_torque_threshold = 150;  % Nm, total decel torque below this gives a meaningless bias
bias_speed_min_mph   = 2;    % mph
bias_accel_max       = -1;   % m/s^2
speed_plot           = linspace(1, 80, 500);  % mph, x-grid for theoretical curves

driveCats   = {'Autocross', 'Endurance', 'Other/Misc'};
driveColors = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19];

%% ================== SELECT DATASETS ==================
[files, path] = uigetfile('*.txt', 'Select driving dataset file(s)', 'MultiSelect', 'on');
if isequal(files, 0)
    error('No files selected.');
end
if ischar(files)
    files = {files};
end
nFiles = numel(files);

%% ================== LOAD, GATE, AND SEGMENT EACH DATASET ==================
datasets = struct('name', {}, 'driveType', {}, 'has_bms', {}, 'summary', {}, 'events', {}, ...
    'bias_brake_samples', {}, 'bias_decel_samples', {}, ...
    'regen_speed_mph', {}, 'regen_Tbias', {}, 'regen_Twheel', {}, 'regen_power_kW', {}, ...
    'regen_bms_current', {}, 'decel_bias_speed_mph', {}, 'decel_bias_value', {});
skipped = {};

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

    if size(raw, 1) < 3
        fprintf('  SKIPPED: fewer than 3 samples.\n');
        skipped{end+1} = sprintf('%s (fewer than 3 samples)', files{k}); %#ok<SAGROW>
        continue
    end

    bms_col = find_bms_column(fname, ncols);
    parsed  = parse_dataset_columns(raw, fmt, bms_col);
    derived = compute_derived_quantities(parsed, VP);
    has_bms = ~isempty(bms_col);

    t  = parsed.t;
    dt = [t(2) - t(1); diff(t)];   % same length as t, first step assumed equal to the second

    % ---- Which samples are regenerating ----
    % With pack current available this is BrakeDataAnalysis.m's own test,
    % unchanged. Without it, the wheel-side regen torque stands in - it
    % measures the same event one step further down the driveline.
    if has_bms
        fprintf('  BMS current: column %d\n', bms_col);
        regen_idx = derived.bms_current < 0 & derived.speed_mph > min_regen_speed;
    else
        fprintf('  BMS current: not present in this file (using regen torque instead)\n');
        regen_idx = derived.total_Twheel_regen < regen_torque_threshold & ...
                    derived.speed_mph > min_regen_speed;
    end

    % ---- Power streams (BrakeDataAnalysis.m physics, unchanged) ----
    decel_idx   = derived.accelx < 0;
    braking_idx = derived.fl_Ttotal < 0 | derived.fr_Ttotal < 0 | ...
                  derived.rl_Ttotal < 0 | derived.rr_Ttotal < 0;

    P_regen = abs(derived.total_regen_power);   % W, already gated to decel
    P_mech  = derived.total_mech_power;         % W, already gated to pressure
    P_aero  = derived.P_aero;                   % W

    % ---- REGEN GATE ----
    % A file with no regen says nothing about the regen contribution, and
    % averaging it in would drag every result toward zero.
    regen_J_gate = sum(P_regen .* dt);
    diss_J_gate  = regen_J_gate + sum(P_mech .* dt) + sum(P_aero .* double(decel_idx) .* dt);
    regen_share  = 100 * regen_J_gate / max(diss_J_gate, eps);
    if has_bms
        keep   = any(derived.bms_current <= min_regen_current_A & derived.speed_mph > min_regen_speed);
        reason = sprintf('pack current never reached %g A above %g mph (min %.1f A)', ...
            min_regen_current_A, min_regen_speed, min(derived.bms_current));
    else
        keep   = regen_share >= min_regen_share_pct;
        reason = sprintf('regen is only %.2f%% of dissipated energy (threshold %g%%)', ...
            regen_share, min_regen_share_pct);
    end
    if ~keep
        fprintf('  SKIPPED: no regen data - %s.\n', reason);
        skipped{end+1} = sprintf('%s (no regen data: %s)', files{k}, reason); %#ok<SAGROW>
        continue
    end
    fprintf('  Regen share of dissipated energy: %.1f%%\n', regen_share);

    driveType = detectDriveType(fname);
    fprintf('  Drive type: %s\n', driveType);

    % ---- Event segmentation ----
    active       = derived.accelx < decel_accel_threshold & derived.speed_mph > min_event_speed_mph;
    event_ranges = find_events(active, t, bridge_gap_s, min_event_dur_s);
    fprintf('  Found %d braking event(s).\n', size(event_ranges, 1));
    if isempty(event_ranges)
        fprintf('  SKIPPED: no braking events detected.\n');
        skipped{end+1} = sprintf('%s (no braking events)', files{k}); %#ok<SAGROW>
        continue
    end

    % ---- Predicted regen at the actual current limit ----
    % Identical model to BrakeDataAnalysis.m: the most restrictive of the
    % current-dependent max-torque curve, the controller power limit, and
    % the absolute torque limit, applied over the regen samples.
    omega_from_speed        = max(derived.speed_mph * 0.44704 / WheelR, 0.1);  % rad/s
    torque_from_power_limit = regen_power_limit ./ omega_from_speed;           % Nm
    max_regen_torque        = regen_max_curve_a * derived.speed_mph.^(regen_max_curve_b) * regen_efficiency;
    max_regen_torque(derived.speed_mph < 1) = 0;
    effective_max_torque    = min(min(max_regen_torque, torque_from_power_limit), torque_limit_total_wheel);
    predicted_regen_torque  = -effective_max_torque .* double(regen_idx);
    predicted_regen_power   = (predicted_regen_torque .* derived.Tbias_regen / 2)     .* derived.fl_omega_wheel + ...
                              (predicted_regen_torque .* derived.Tbias_regen / 2)     .* derived.fr_omega_wheel + ...
                              (predicted_regen_torque .* (1-derived.Tbias_regen) / 2) .* derived.rl_omega_wheel + ...
                              (predicted_regen_torque .* (1-derived.Tbias_regen) / 2) .* derived.rr_omega_wheel;

    % ---- Per-dataset energy totals (J) ----
    S = struct();
    S.regen_J      = sum(P_regen .* dt);
    S.mech_J       = sum(P_mech  .* dt);
    S.aero_decel_J = sum(P_aero .* double(decel_idx) .* dt);
    S.aero_regen_J = sum(P_aero .* double(regen_idx) .* dt);
    S.pred_regen_J = sum(abs(predicted_regen_power .* dt));
    S.aero_open_J       = sum(derived.F_aero_open   .* max(derived.velx, 0) .* double(derived.drs_state == 1) .* dt);
    S.aero_closed_J     = sum(derived.F_aero_closed .* max(derived.velx, 0) .* double(derived.drs_state == 0) .* dt);
    S.aero_all_open_J   = sum(derived.F_aero_open   .* max(derived.velx, 0) .* dt);
    S.aero_all_closed_J = sum(derived.F_aero_closed .* max(derived.velx, 0) .* dt);
    S.n_events          = size(event_ranges, 1);

    % ---- Per-event records ----
    events = struct('dataset', {}, 'driveType', {}, 'duration_s', {}, ...
        'avg_decel_g', {}, 'peak_decel_g', {}, 'start_speed_mph', {}, 'end_speed_mph', {}, ...
        'E_regen_kJ', {}, 'E_fric_kJ', {}, 'E_aero_kJ', {}, 'E_total_kJ', {}, 'regen_frac_pct', {});

    for e = 1:size(event_ranges, 1)
        i0 = event_ranges(e, 1);
        i1 = event_ranges(e, 2);
        duration_s = t(i1) - t(i0);
        if duration_s <= 0
            continue
        end
        idx = i0:i1;

        E_regen_kJ = sum(P_regen(idx) .* dt(idx)) / 1000;
        E_fric_kJ  = sum(P_mech(idx)  .* dt(idx)) / 1000;
        E_aero_kJ  = sum(P_aero(idx)  .* dt(idx)) / 1000;
        E_total_kJ = E_regen_kJ + E_fric_kJ + E_aero_kJ;

        events(end+1) = struct('dataset', files{k}, 'driveType', driveType, ...
            'duration_s', duration_s, ...
            'avg_decel_g',  (derived.velx(i0) - derived.velx(i1)) / duration_s / 9.81, ...
            'peak_decel_g', abs(min(derived.accelx(idx))) / 9.81, ...
            'start_speed_mph', derived.speed_mph(i0), 'end_speed_mph', derived.speed_mph(i1), ...
            'E_regen_kJ', E_regen_kJ, 'E_fric_kJ', E_fric_kJ, 'E_aero_kJ', E_aero_kJ, ...
            'E_total_kJ', E_total_kJ, ...
            'regen_frac_pct', 100 * E_regen_kJ / max(E_total_kJ, eps)); %#ok<SAGROW>
    end

    % ---- Sample pools for the figures ----
    in_event = false(size(t));
    for e = 1:size(event_ranges, 1)
        in_event(event_ranges(e,1):event_ranges(e,2)) = true;
    end

    % Figure 1 uses only samples where the friction brakes are actually
    % applied - otherwise Tbias_brake is an identically-zero placeholder
    % (0/0 guarded to 0) and would pile a false spike at x = 0.
    bias_brake_idx = in_event & derived.mech_braking_idx;
    % Figure 2 likewise needs a meaningful denominator: below
    % min_torque_threshold the front/total ratio is noise divided by noise.
    bias_decel_idx = in_event & braking_idx & ...
                     abs(derived.ftot_Ttotal + derived.rtot_Ttotal) > min_torque_threshold;

    % Figure 4 keeps BrakeDataAnalysis.m's original stricter gating.
    decel_bias_plot_idx = braking_idx & ...
                          abs(derived.ftot_Ttotal + derived.rtot_Ttotal) > min_torque_threshold & ...
                          derived.speed_mph > bias_speed_min_mph & ...
                          derived.accelx < bias_accel_max;

    if has_bms
        regen_bms_current = derived.bms_current(regen_idx);
    else
        regen_bms_current = [];
    end

    datasets(end+1) = struct('name', files{k}, 'driveType', driveType, 'has_bms', has_bms, ...
        'summary', S, 'events', events, ...
        'bias_brake_samples', derived.Tbias_brake(bias_brake_idx), ...
        'bias_decel_samples', derived.Tbias_decel(bias_decel_idx), ...
        'regen_speed_mph',   derived.speed_mph(regen_idx), ...
        'regen_Tbias',       derived.Tbias_regen(regen_idx), ...
        'regen_Twheel',      abs(derived.total_Twheel_regen(regen_idx)), ...
        'regen_power_kW',    P_regen(regen_idx) / 1000, ...
        'regen_bms_current', regen_bms_current, ...
        'decel_bias_speed_mph', derived.speed_mph(decel_bias_plot_idx), ...
        'decel_bias_value',     derived.Tbias_decel(decel_bias_plot_idx)); %#ok<SAGROW>
end

if isempty(datasets)
    error('No usable datasets loaded (every selected file was skipped).');
end

fprintf('\n%d of %d selected dataset(s) used.\n', numel(datasets), nFiles);
if ~isempty(skipped)
    fprintf('Excluded:\n');
    for s = 1:numel(skipped)
        fprintf('  - %s\n', skipped{s});
    end
end

allEvents = [datasets.events];
allTypes  = {allEvents.driveType};

%% ================== CONSOLE OUTPUT ==================
fprintf('\n\n################ PER-DATASET ################\n');
for k = 1:numel(datasets)
    print_energy_summary(sprintf('%s  [%s]', datasets(k).name, datasets(k).driveType), ...
        datasets(k).summary, current_limit_actual, min_energy_threshold);
end

fprintf('\n\n################ PER-DRIVE-TYPE ################\n');
for c = 1:numel(driveCats)
    memberIdx = strcmp({datasets.driveType}, driveCats{c});
    if ~any(memberIdx)
        fprintf('\n===== %s: no datasets =====\n', driveCats{c});
        continue
    end
    print_energy_summary(sprintf('%s  (%d dataset(s))', driveCats{c}, sum(memberIdx)), ...
        sum_summaries([datasets(memberIdx).summary]), current_limit_actual, min_energy_threshold);
end

fprintf('\n\n################ ALL DATASETS POOLED ################\n');
totalSummary = sum_summaries([datasets.summary]);
print_energy_summary(sprintf('ALL (%d datasets)', numel(datasets)), ...
    totalSummary, current_limit_actual, min_energy_threshold);

%% ---- HEADLINE: REGEN CONTRIBUTION ----
% Two different averages, both worth having. The energy-weighted figure is
% total regen over total dissipated - it answers "what share of the energy
% came back". The mean-per-event figure averages each event's own share
% equally - it answers "in a typical braking event, how much comes back".
% Big stops dominate the first; nothing dominates the second.
fprintf('\n\n################ REGEN CONTRIBUTION ################\n');
fprintf('Denominator is total dissipated energy: regen + friction brakes + aero drag.\n\n');
fprintf('%-16s %18s %28s\n', 'Group', 'Energy-weighted', 'Mean per event');
fprintf('%s\n', repmat('-', 1, 66));
for c = 1:numel(driveCats)
    memberIdx = strcmp({datasets.driveType}, driveCats{c});
    if ~any(memberIdx)
        continue
    end
    print_regen_contribution(driveCats{c}, sum_summaries([datasets(memberIdx).summary]), ...
        allEvents(strcmp(allTypes, driveCats{c})));
end
fprintf('%s\n', repmat('-', 1, 66));
print_regen_contribution('ALL', totalSummary, allEvents);

%% ================================================================
%  FIGURES
% ================================================================

%% ---- Figure 1: Mechanical brake torque bias distribution ----
% Every in-event sample pooled across all datasets. The axis is clipped to
% [0 1] for display only - no sample is dropped from the fit or the stats.
figure('Name', 'Mechanical Brake Torque Bias Distribution');
biasBrake = vertcat(datasets.bias_brake_samples);
fitBrake  = plot_distribution_overlay({biasBrake}, {'Brake Torque Bias'}, [0 0 1], ...
    'Mechanical Brake Torque Bias (Front/Total)', ...
    sprintf('Mechanical Brake Torque Bias during Braking Events (n = %d samples)', numel(biasBrake)));
xlim([0 1]);

%% ---- Figure 2: Total decel torque bias distribution ----
figure('Name', 'Total Decel Torque Bias Distribution');
biasDecel = vertcat(datasets.bias_decel_samples);
fitDecel  = plot_distribution_overlay({biasDecel}, {'Total Decel Torque Bias'}, [0 0 0], ...
    'Total Decel Torque Bias (Front/Total)', ...
    sprintf('Total Decel Torque Bias during Braking Events (n = %d samples)', numel(biasDecel)));
xlim([0 1]);

fprintf('\n\n################ TORQUE BIAS DISTRIBUTIONS ################\n');
report_distribution('Mechanical brake torque bias', fitBrake{1});
report_distribution('Total decel torque bias',      fitDecel{1});

%% ---- Figure 3: Regen torque bias vs speed ----
figure('Name', 'Regen Torque Bias vs Speed');
scatter(vertcat(datasets.regen_speed_mph), vertcat(datasets.regen_Tbias), 50, 'k.');
xlabel('Vehicle Speed (mph)'); ylabel('Regen Torque Bias (Front/Total)');
title('Regen Torque Bias vs Vehicle Speed'); grid on;

%% ---- Figure 4: Total decel torque bias vs speed ----
figure('Name', 'Total Decel Torque Bias vs Speed');
scatter(vertcat(datasets.decel_bias_speed_mph), vertcat(datasets.decel_bias_value), 50, 'k.');
ylim([0 1]);
xlabel('Vehicle Speed (mph)'); ylabel('Total Decel Torque Bias (Front/Total)');
title('Total Decel Torque Bias vs Vehicle Speed during Braking Events'); grid on;

%% ---- Figure 5: Regen total wheel torque vs speed ----
% Split in two because the color axis is not available for every file.
% Datasets WITH pack current keep BrakeDataAnalysis.m's original coloring;
% datasets without it are colored by regen power, the closest available
% stand-in and the quantity the controller power limit acts on. Mixing the
% two into one colorbar would put amps and kilowatts on the same scale.
regenSpeed  = vertcat(datasets.regen_speed_mph);
regenTwheel = vertcat(datasets.regen_Twheel);
bmsMask     = [datasets.has_bms];

if any(bmsMask)
    figure('Name', 'Regen Total Wheel Torque vs Speed (BMS current)');
    scatter(vertcat(datasets(bmsMask).regen_speed_mph), ...
            vertcat(datasets(bmsMask).regen_Twheel), 20, ...
            vertcat(datasets(bmsMask).regen_bms_current), '.');
    colormap('jet');
    clabel = colorbar; clabel.Label.String = 'Current (A)';
    xlabel('Vehicle Speed (mph)'); ylabel('Regen Wheel Torque (Nm)');
    title(sprintf('Regen Total Wheel Torque vs Speed (%d dataset(s) with BMS current)', ...
        sum(bmsMask)));
    grid on;
end

if any(~bmsMask)
    figure('Name', 'Regen Total Wheel Torque vs Speed (regen power)');
    scatter(vertcat(datasets(~bmsMask).regen_speed_mph), ...
            vertcat(datasets(~bmsMask).regen_Twheel), 20, ...
            vertcat(datasets(~bmsMask).regen_power_kW), '.');
    colormap('jet');
    clabel = colorbar; clabel.Label.String = 'Regen Power (kW)';
    xlabel('Vehicle Speed (mph)'); ylabel('Regen Wheel Torque (Nm)');
    title(sprintf('Regen Total Wheel Torque vs Speed (%d dataset(s) without BMS current)', ...
        sum(~bmsMask)));
    grid on;
end

%% ---- Figure 6: Theoretical max regen torque vs speed at current limits ----
figure('Name', 'Theoretical Regen Torque vs Speed');
omega_plot        = max(speed_plot * 0.44704 / WheelR, 0.1);   % rad/s from mph
torque_power_plot = regen_power_limit ./ omega_plot;           % power-limit curve, Nm

for i = 1:numel(current_limits)
    scale_factor    = current_limits(i) / current_limit_actual;
    max_regen_speed = regen_max_curve_a * speed_plot.^(regen_max_curve_b) * scale_factor;
    effective_limit = min(min(max_regen_speed, torque_power_plot), torque_limit_total_wheel);
    plot(speed_plot, effective_limit, limit_colors{mod(i-1, numel(limit_colors)) + 1}, ...
        'LineWidth', 2, 'DisplayName', sprintf('%dA limit', current_limits(i))); hold on;
end
plot(speed_plot, torque_power_plot, 'c--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Power Limit (%.0f kW)', regen_power_limit/1000));
scatter(regenSpeed, regenTwheel, 8, [0.5 0.5 0.5], '.', 'DisplayName', 'Measured Regen Torque');
yline(torque_limit_total_wheel, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Torque Limit');
xlabel('Vehicle Speed (mph)'); ylabel('Max Regen Torque (Nm)');
title('Theoretical Max Regen Torque vs Speed at Different Current Limits');
legend('Location', 'best'); grid on;
ylim([0, torque_limit_total_wheel * 1.2]);

%% ---- Figure 7: Cumulative braking energy breakdown by drive type ----
figure('Name', 'Cumulative Braking Energy Breakdown by Drive Type');
for c = 1:numel(driveCats)
    subplot(1, 3, c);
    memberIdx = strcmp({datasets.driveType}, driveCats{c});
    if ~any(memberIdx)
        axis off;
        title(sprintf('%s\n(no datasets)', driveCats{c}));
        continue
    end
    Sc   = sum_summaries([datasets(memberIdx).summary]);
    vals = [Sc.mech_J, Sc.regen_J, Sc.aero_decel_J] / 1000;   % kJ
    names = {'Friction', 'Regen', 'Aero'};
    % The kJ/percent detail goes in the legend rather than in wedge labels:
    % as text labels they run into the title and into the neighbouring pie.
    pie(vals);
    legend(arrayfun(@(v, n) sprintf('%s: %.0f kJ (%.1f%%)', n{1}, v, 100*v/max(sum(vals), eps)), ...
        vals, names, 'UniformOutput', false), ...
        'Location', 'southoutside', 'FontSize', 7);
    title(sprintf('%s\n%d events, %.0f kJ total', driveCats{c}, ...
        sum(strcmp(allTypes, driveCats{c})), sum(vals)));
end
sgtitle('Cumulative Braking Energy Breakdown');

%% ---- Figure 8: Per-event energy distributions by drive type ----
figure('Name', 'Per-Event Braking Energy Distributions by Drive Type');
energyColors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];   % regen, friction, aero
fprintf('\n\n################ PER-EVENT ENERGY DISTRIBUTIONS ################\n');
for c = 1:numel(driveCats)
    subplot(1, 3, c);
    mask = strcmp(allTypes, driveCats{c});
    if ~any(mask)
        axis off;
        title(sprintf('%s\n(no datasets)', driveCats{c}));
        continue
    end
    ev   = allEvents(mask);
    fits = plot_distribution_overlay({[ev.E_regen_kJ], [ev.E_fric_kJ], [ev.E_aero_kJ]}, ...
        {'Regen', 'Friction', 'Aero'}, energyColors, 'Energy per Braking Event (kJ)', ...
        sprintf('%s (n = %d events)', driveCats{c}, numel(ev)));
    % Display-only x clip. A few very large stops otherwise stretch the axis
    % far enough that the bulk of the events collapses against the origin.
    % Nothing is excluded from the histogram, the fits, or the stats.
    xHigh = prctile([ev.E_regen_kJ, ev.E_fric_kJ, ev.E_aero_kJ], 99);
    if isfinite(xHigh) && xHigh > 0
        xlim([0, xHigh]);
    end
    fprintf('\n--- %s (n = %d events) ---\n', driveCats{c}, numel(ev));
    report_distribution('Regen',    fits{1});
    report_distribution('Friction', fits{2});
    report_distribution('Aero',     fits{3});
end
sgtitle('Per-Event Energy Distributions (regen / friction / aero)');

%% ---- Figure 9: Regen energy vs decel rate and event duration ----
figure('Name', 'Regen Energy Correlations');
corrSpecs = {'avg_decel_g', 'Average Deceleration (g)'; ...
             'duration_s',  'Braking Event Duration (s)'};
for p = 1:size(corrSpecs, 1)
    subplot(1, 2, p); hold on; grid on;
    for c = 1:numel(driveCats)
        mask = strcmp(allTypes, driveCats{c});
        if ~any(mask)
            continue
        end
        ev = allEvents(mask);
        scatter([ev.(corrSpecs{p,1})], [ev.E_regen_kJ], 20, driveColors(c,:), 'filled', ...
            'MarkerFaceAlpha', 0.55, 'DisplayName', driveCats{c});
    end
    xlabel(corrSpecs{p,2}); ylabel('Regen Energy per Event (kJ)');
    title(sprintf('Regen Energy vs %s', corrSpecs{p,2}), 'Interpreter', 'none');
    if p == 1
        legend('Location', 'best');
    end
end
sgtitle('Per-Event Regen Energy Correlations');


%% ================================================================
%  LOCAL FUNCTIONS
% ================================================================

function bms_col = find_bms_column(fname, ncols)
% Locates the pack-current channel BY NAME in the channel-name header row,
% returning [] when the file has no such channel.
%
% Column position is not reliable here and column count does not imply a
% layout: a 22-column RAW logger file carries BMS current in column 7,
% while a 22-column PREPROCESSED export carries Steering_pct there. The
% channels the positional parser actually reads (pressures, motor torque/
% velocity, velocity, accel) sit at the same indices in both, so only pack
% current needs finding by name.
bms_col = [];
fid = fopen(fname, 'r');
if fid < 0
    return
end
cleanup = onCleanup(@() fclose(fid));
for lineNum = 1:6
    line = fgetl(fid);
    if ~ischar(line)
        return
    end
    names = strsplit(line, sprintf('\t'));
    if numel(names) ~= ncols
        continue   % not the channel-name row
    end
    hit = find(contains(lower(names), 'bms'), 1);
    if ~isempty(hit)
        bms_col = hit;
        return
    end
end
end


function parsed = parse_dataset_columns(data, fmt, bms_col)
% Local copy of padSpecificPower.m's parser, EXTENDED with accelx (the
% event detector) and with pack current, which is read from the column
% find_bms_column located by name rather than by a fixed position.
parsed.t           = data(:,1);
parsed.fl_temp_adc = data(:,2);
parsed.fr_temp_adc = data(:,3);
parsed.rl_temp_adc = data(:,4);
parsed.rr_temp_adc = data(:,5);

if isempty(bms_col)
    parsed.bms_current = [];
else
    parsed.bms_current = data(:, bms_col);
end

switch fmt
    case {'drs_imu', 'drs'}
        parsed.drs_state         = data(:,6);
        parsed.frontpressure_adc = data(:,8);
        parsed.rearpressure_adc  = data(:,12);
        parsed.fl_Tmotor_Mn = data(:,13); parsed.fl_vmotor = data(:,14);
        parsed.fr_Tmotor_Mn = data(:,15); parsed.fr_vmotor = data(:,16);
        parsed.rl_Tmotor_Mn = data(:,17); parsed.rl_vmotor = data(:,18);
        parsed.rr_Tmotor_Mn = data(:,19); parsed.rr_vmotor = data(:,20);
        if strcmp(fmt, 'drs_imu')
            parsed.velx   = data(:,26);
            parsed.accelx = data(:,27);
        else
            parsed.velx   = data(:,21);
            parsed.accelx = data(:,22);
        end
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
% Local copy of padSpecificPower.m's compute_derived_quantities (same
% per-single-component basis convention - see loadVehicleParams.m's header),
% EXTENDED with the wheel/total/regen torque quantities and the mechanical
% braking power stream that BrakeDataAnalysis.m computes.
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

fl_Twheel = fl_Tmotor * VP.gear_ratio;  fr_Twheel = fr_Tmotor * VP.gear_ratio;
rl_Twheel = rl_Tmotor * VP.gear_ratio;  rr_Twheel = rr_Tmotor * VP.gear_ratio;

fl_vwheel = parsed.fl_vmotor / VP.gear_ratio;  fr_vwheel = parsed.fr_vmotor / VP.gear_ratio;
rl_vwheel = parsed.rl_vmotor / VP.gear_ratio;  rr_vwheel = parsed.rr_vmotor / VP.gear_ratio;

fl_omega_wheel = fl_vwheel * (2*pi/60);  fr_omega_wheel = fr_vwheel * (2*pi/60);
rl_omega_wheel = rl_vwheel * (2*pi/60);  rr_omega_wheel = rr_vwheel * (2*pi/60);
fl_omega_motor = parsed.fl_vmotor * (2*pi/60);  fr_omega_motor = parsed.fr_vmotor * (2*pi/60);
rl_omega_motor = parsed.rl_vmotor * (2*pi/60);  rr_omega_motor = parsed.rr_vmotor * (2*pi/60);

% Clamp force uses the SINGLE-SIDE piston area and the MEAN PAD RADIUS as
% the friction lever arm - see loadVehicleParams.m for why.
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

% Total decel torque = motor (regen) + friction brake, per corner
fl_Ttotal = fl_Twheel + fl_Tbrake;  fr_Ttotal = fr_Twheel + fr_Tbrake;
rl_Ttotal = rl_Twheel + rl_Tbrake;  rr_Ttotal = rr_Twheel + rr_Tbrake;
ftot_Ttotal = fl_Ttotal + fr_Ttotal;
rtot_Ttotal = rl_Ttotal + rr_Ttotal;
Tbias_decel = ftot_Ttotal ./ (ftot_Ttotal + rtot_Ttotal);
Tbias_decel(isnan(Tbias_decel)) = 0;

% Regen-only wheel torque and its front/rear split
ftot_Twheel_regen = min(fl_Twheel, 0) + min(fr_Twheel, 0);
rtot_Twheel_regen = min(rl_Twheel, 0) + min(rr_Twheel, 0);
Tbias_regen = ftot_Twheel_regen ./ (ftot_Twheel_regen + rtot_Twheel_regen);
Tbias_regen(isnan(Tbias_regen)) = 0;

F_aero_open   = VP.aero_open_a   * speed_mph.^2 + VP.aero_open_b   * speed_mph + VP.aero_open_c;
F_aero_closed = VP.aero_closed_a * speed_mph.^2 + VP.aero_closed_b * speed_mph + VP.aero_closed_c;
F_aero = F_aero_open .* double(parsed.drs_state == 1) + F_aero_closed .* double(parsed.drs_state == 0);
F_aero = max(F_aero, 0);

% Regen power, per corner, only while decelerating
decelerating_idx = parsed.accelx < 0;
fl_regen_power = min(fl_Tmotor .* fl_omega_motor, 0) .* double(decelerating_idx);
fr_regen_power = min(fr_Tmotor .* fr_omega_motor, 0) .* double(decelerating_idx);
rl_regen_power = min(rl_Tmotor .* rl_omega_motor, 0) .* double(decelerating_idx);
rr_regen_power = min(rr_Tmotor .* rr_omega_motor, 0) .* double(decelerating_idx);

% Mechanical braking power, only while the pads are actually engaged
mech_braking_idx = (frontpressure > VP.min_pressure | rearpressure > VP.min_pressure);
total_mech_power = (abs(fl_Tbrake .* fl_omega_wheel) + abs(fr_Tbrake .* fr_omega_wheel) + ...
                    abs(rl_Tbrake .* rl_omega_wheel) + abs(rr_Tbrake .* rr_omega_wheel)) ...
                   .* double(mech_braking_idx);

derived.velx          = velx;
derived.speed_mph     = speed_mph;
derived.accelx        = parsed.accelx;
derived.bms_current   = parsed.bms_current;
derived.drs_state     = parsed.drs_state;
derived.frontpressure = frontpressure;
derived.rearpressure  = rearpressure;
derived.fr_temp_F     = fr_temp_F;
derived.rr_temp_F     = rr_temp_F;

derived.fl_omega_wheel = fl_omega_wheel;  derived.fr_omega_wheel = fr_omega_wheel;
derived.rl_omega_wheel = rl_omega_wheel;  derived.rr_omega_wheel = rr_omega_wheel;

derived.fl_Ttotal = fl_Ttotal;  derived.fr_Ttotal = fr_Ttotal;
derived.rl_Ttotal = rl_Ttotal;  derived.rr_Ttotal = rr_Ttotal;
derived.ftot_Ttotal = ftot_Ttotal;  derived.rtot_Ttotal = rtot_Ttotal;

derived.Tbias_brake = Tbias_brake;
derived.Tbias_decel = Tbias_decel;
derived.Tbias_regen = Tbias_regen;
derived.total_Twheel_regen = ftot_Twheel_regen + rtot_Twheel_regen;

derived.F_aero        = F_aero;
derived.F_aero_open   = F_aero_open;
derived.F_aero_closed = F_aero_closed;
derived.P_aero        = F_aero .* velx;

derived.front_axle_regen_power = fl_regen_power + fr_regen_power;
derived.rear_axle_regen_power  = rl_regen_power + rr_regen_power;
derived.total_regen_power      = derived.front_axle_regen_power + derived.rear_axle_regen_power;

derived.mech_braking_idx = mech_braking_idx;
derived.total_mech_power = total_mech_power;
end


function event_ranges = find_events(active, t, bridge_gap_s, min_event_dur_s)
% Local copy of padSpecificPower.m's find_events. Run-length-encodes the
% `active` logical vector into [start_idx end_idx] ranges, bridges gaps
% shorter than bridge_gap_s, then drops events shorter than min_event_dur_s.

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
% Local copy of padSpecificPower.m's detectDriveType.
%
% Primary signal: the FILENAME itself. Checked first because RawFormatOnly
% exports from F34BrakeDataPreprocessor.m always write "Segment type:
% raw_format_only" in their header regardless of the file's true drive type.
[~, baseName, ~] = fileparts(fname);
lowerName = lower(baseName);
if contains(lowerName, 'endurance')
    driveType = 'Endurance';
    return
end
if contains(lowerName, 'autocross') || contains(lowerName, 'autox')
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
            if contains(kind, 'endurance')
                driveType = 'Endurance';
            elseif contains(kind, 'autocross')
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


function S = sum_summaries(Sarr)
% Pools an array of per-dataset energy summaries by summing every field.
% Every field is an extensive quantity (joules or a count), so summing is
% the correct pooling operation for all of them - percentages and ratios
% are always recomputed from the pooled totals, never averaged.
S = Sarr(1);
flds = fieldnames(S);
for f = 1:numel(flds)
    S.(flds{f}) = sum([Sarr.(flds{f})]);
end
end


function print_energy_summary(label, S, current_limit_actual, min_energy_threshold)
% The BrakeDataAnalysis.m console blocks, minus the temperature output.
total_braking_J = S.mech_J + S.regen_J;
total_decel_J   = total_braking_J + S.aero_decel_J;

fprintf('\n===== %s =====\n', label);

fprintf('\n--- Braking Energy Summary ---\n');
fprintf('Total Braking Energy:   %.2f kJ  (%.2f Wh)\n', total_braking_J/1000, total_braking_J/3600);
fprintf('Regenerative Braking:   %.2f kJ  (%.2f Wh)\n', S.regen_J/1000, S.regen_J/3600);
fprintf('Mechanical Braking:     %.2f kJ  (%.2f Wh)\n', S.mech_J/1000,  S.mech_J/3600);
fprintf('Regen Fraction:         %.1f%%\n', 100 * S.regen_J / max(total_braking_J, eps));
fprintf('Number of Braking Events: %d\n', S.n_events);

fprintf('\n--- Aero Drag During Regen ---\n');
fprintf('Aero Energy during regen: %.2f kJ\n', S.aero_regen_J/1000);

fprintf('\n--- %dA Regen Model vs Actual ---\n', current_limit_actual);
fprintf('Actual Regen Energy:    %.2f kJ\n', S.regen_J/1000);
fprintf('Predicted Regen Energy: %.2f kJ\n', S.pred_regen_J/1000);
if S.regen_J/1000 > min_energy_threshold
    fprintf('Percent Difference:     %.2f%%\n', 100 * (S.pred_regen_J - S.regen_J) / S.regen_J);
else
    fprintf('Percent Difference:     n/a (actual regen energy below the %.0f kJ threshold)\n', ...
        min_energy_threshold);
end
fprintf('Aero fraction of predicted: %.1f%%\n', 100 * S.aero_regen_J / max(S.pred_regen_J, eps));

fprintf('\n--- Braking Energy Breakdown ---\n');
fprintf('Mechanical:   %.2f kJ  (%.1f%%)\n', S.mech_J/1000,       100 * S.mech_J       / max(total_decel_J, eps));
fprintf('Regenerative: %.2f kJ  (%.1f%%)\n', S.regen_J/1000,      100 * S.regen_J      / max(total_decel_J, eps));
fprintf('Aero:         %.2f kJ  (%.1f%%)\n', S.aero_decel_J/1000, 100 * S.aero_decel_J / max(total_decel_J, eps));
fprintf('Total:        %.2f kJ\n', total_decel_J/1000);

fprintf('\n--- DRS Drag Energy Analysis ---\n');
actual_aero_J = S.aero_open_J + S.aero_closed_J;
fprintf('Actual aero energy (DRS open):      %.2f kJ\n', S.aero_open_J/1000);
fprintf('Actual aero energy (DRS closed):    %.2f kJ\n', S.aero_closed_J/1000);
fprintf('Total actual aero energy:           %.2f kJ\n', actual_aero_J/1000);
fprintf('If always DRS open:                 %.2f kJ\n', S.aero_all_open_J/1000);
fprintf('If always DRS closed:               %.2f kJ\n', S.aero_all_closed_J/1000);
fprintf('DRS energy saving vs always closed: %.2f kJ (%.1f%%)\n', ...
    (S.aero_all_closed_J - actual_aero_J)/1000, ...
    100 * (S.aero_all_closed_J - actual_aero_J) / max(S.aero_all_closed_J, eps));
end


function print_regen_contribution(label, S, events)
% The headline answer, both ways: energy-weighted (big stops dominate) and
% averaged equally over events (a typical stop).
total_decel_J  = S.mech_J + S.regen_J + S.aero_decel_J;
energyWeighted = 100 * S.regen_J / max(total_decel_J, eps);

fracs = [events.regen_frac_pct];
if isempty(fracs)
    fprintf('%-16s %17.1f%% %28s\n', label, energyWeighted, 'no events');
else
    fprintf('%-16s %17.1f%% %16.1f%% +/- %4.1f%%  (n = %d)\n', ...
        label, energyWeighted, mean(fracs), std(fracs), numel(fracs));
end
end


function report_distribution(label, fitResult)
% One-line console summary of a fitted distribution.
if isempty(fitResult.bestParams)
    paramsStr = '-';
else
    paramsStr = sprintf('%.4g, ', fitResult.bestParams);
    paramsStr = paramsStr(1:end-2);   % drop trailing ", "
end
fprintf('%-28s n=%-7d mean=%-9.3f median=%-9.3f std=%-9.3f skew=%-8.3f best fit: %s (%s)\n', ...
    label, fitResult.n, fitResult.mean, fitResult.median, fitResult.std, fitResult.skew, ...
    fitResult.bestName, paramsStr);
end


function fitResult = fit_distribution_summary(data)
% Local copy of padSpecificPower.m's fit_distribution_summary. Fits
% candidate distributions to `data` and selects the best by AIC
% (2*k - 2*loglikelihood), computed manually since fitdist objects don't
% expose a uniform AIC property across distribution types. Requires the
% Statistics and Machine Learning Toolbox.
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

% Exponential is included alongside padSpecificPower.m's four candidates
% because per-event energies are zero-inflated and heavily right-skewed - an
% event braked purely on regen contributes exactly 0 friction energy. The
% strictly-positive candidates (Lognormal, Gamma, Weibull) degrade or fail
% outright on those zeros, which would otherwise hand every such series to
% Normal by default despite an obviously non-normal shape.
candidates = {'Normal', 'Lognormal', 'Gamma', 'Weibull', 'Exponential'};
bestAIC = Inf;

% A candidate that cannot represent the data is handled by losing on AIC
% below; its fitting warnings add nothing and would bury the real output.
warnState  = warning('off', 'all');
restoreWarn = onCleanup(@() warning(warnState));

for c = 1:numel(candidates)
    try
        pd = fitdist(data, candidates{c});
        % Clip pdf values away from exact 0 before taking log so a single
        % far-outlier point penalizes a fit instead of disqualifying it.
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


function fits = plot_distribution_overlay(dataSets, labels, colorList, xLabelStr, titleStr)
% Generalization of padSpecificPower.m's front/rear plot_distribution_overlay
% to an arbitrary number of labeled series (Figure 8 needs three). Draws
% probability-density histograms with each series' best-fit PDF on top, into
% the CURRENT axes. Returns each series' fit so the caller can report the
% stats without refitting.
fits = cell(1, numel(dataSets));
hold on; grid on;

maxDensity = 0;
for s = 1:numel(dataSets)
    d = dataSets{s}(:);
    d = d(isfinite(d));
    fits{s} = fit_distribution_summary(d);
    if isempty(d)
        continue
    end
    h = histogram(d, 'Normalization', 'pdf', 'FaceColor', colorList(s,:), 'FaceAlpha', 0.35, ...
        'EdgeColor', 'none', 'DisplayName', labels{s});
    maxDensity = max(maxDensity, max(h.Values));
    if ~isempty(fits{s}.pd)
        % The fit curve is drawn in its series' own color and needs no
        % legend row of its own - naming the fit on the histogram's row
        % keeps a 3-series legend from swamping the axes.
        plot(xg_for(d), pdf(fits{s}.pd, xg_for(d)), '-', 'Color', colorList(s,:), ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
        h.DisplayName = sprintf('%s (%s fit)', labels{s}, fits{s}.bestName);
    end
end

% Scale to the histogram, not to the fitted curves: a fit with a pole at
% the left edge (common for zero-inflated per-event energies) would
% otherwise squash the actual data into the bottom of the axes.
if maxDensity > 0
    ylim([0, 1.25 * maxDensity]);
end

xlabel(xLabelStr); ylabel('Probability Density');
title(titleStr, 'Interpreter', 'none');
legend('Location', 'northeast', 'FontSize', 7);
end


function xg = xg_for(d)
% Evaluation grid for a fitted PDF over the span of its data.
xg = linspace(min(d), max(d), 200);
end
