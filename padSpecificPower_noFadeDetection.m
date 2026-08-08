%% ================================================================
%  Brake Pad Specific Power (Heat Flux) Analysis
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
%  Outputs:
%    1) Specific power (W/cm^2) vs. time-since-event-start, for every
%       braking event found, front and rear pads plotted separately.
%    2) The average specific power across all events (front & rear).
%    3) A box-and-whisker plot of each event's average specific power
%       (front vs. rear).
%
%  Datasets are selected via a multi-select file picker, same
%  21/22/27-column auto-detected formats as BrakeDataAnalysis.m.
% ================================================================

clc; clear; close all

%% ================== CALIBRATED MODEL CONSTANTS ==================
% >>> Replace these with your finalized fit from BrakeCoeffOptimizer.m <<<

% Rotor cooling coefficient: h_w(v) = x1*v + b1   [W/m^2-K]
x1_f = 0.055860;
b1_f = 32.025052;

x1_r = 2.682170;
b1_r = 28.826248;

% Pad energy fraction: PadFrac = f(T_rotor [K], P_applied [psi])
% Default below matches the "linear in T only" baseline model.
% If you settled on a pressure-dependent model instead, edit padfrac_fun
% and the coefficients to match (see commented examples).
p(1) = -0.002440977;
p(2) = 2.8126283e-06;
p(3) = -0.00068110755;
p(4) = 1.2826624;
padfrac_fun = @(T, P) p(1).*T + p(2).*T.^2 + p(3).*P + p(4);

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
events_front = struct('dataset', {}, 't_rel', {}, 'q_inst', {}, 'avg_q', {}, 't_decel', {});
events_rear  = struct('dataset', {}, 't_rel', {}, 'q_inst', {}, 'avg_q', {}, 't_decel', {});

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
        x1_f, b1_f, padfrac_fun, derived.total_regen_power, Edrag, ...
        derived.fl_omega_wheel, derived.fr_omega_wheel, ...
        VehicleMass, RotorMass_front, RotorArea_front, I, WheelR, TambK, ...
        A_pad_front_cm2, BrakeFrac, CalibrationFactor, min_pressure, 0.5, TambC);

    sim_rear = simulate_pad_power(t, derived.velx, derived.rearpressure, 1 - derived.Tbias_brake, ...
        x1_r, b1_r, padfrac_fun, derived.total_regen_power, Edrag, ...
        derived.rl_omega_wheel, derived.rr_omega_wheel, ...
        VehicleMass, RotorMass_rear, RotorArea_rear, I, WheelR, TambK, ...
        A_pad_rear_cm2, BrakeFrac, CalibrationFactor, min_pressure, 0.5, TambC);

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
            'q_inst', q_inst_f, 'avg_q', avg_q_f, 't_decel', t_decel); %#ok<SAGROW>
        events_rear(end+1) = struct('dataset', files{k}, 't_rel', t_rel_r, ...
            'q_inst', q_inst_r, 'avg_q', avg_q_r, 't_decel', t_decel); %#ok<SAGROW>
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
    min_pressure, corner_split, TambC) %#ok<INUSD>
% Walks the full timeseries computing rotor temp (needed for PadFrac's
% temperature dependence), per-step pad energy, and instantaneous
% specific power (W/cm^2). PadFrac uses the SAME (1-PadFrac)/PadFrac
% split convention as brake_temp_sim.m: PadFrac is the pad's own share,
% (1-PadFrac) is the rotor's share.

n = length(t);
RotorTempArrayK = zeros(n, 1);
RotorTempArrayK(1) = TambC + 273.15;   % assume ambient at dataset start
pad_energy = zeros(n, 1);   % J, per single pad, per step
q_inst     = zeros(n, 1);   % W/cm^2, per single pad, per step

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
    else
        qout = (prevTemp - TambK) / Rrotor;
        Eout = qout * tbrake;
        deltaTKout = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = prevTemp - deltaTKout;
        % pad_energy(i) and q_inst(i) remain 0 - no engagement this step
    end
end

sim.RotorTempArrayK = RotorTempArrayK;
sim.pad_energy      = pad_energy;
sim.q_inst          = q_inst;
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
