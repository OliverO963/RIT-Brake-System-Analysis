%% F34 Brakes Data Analysis
clc; clear; close all

%% ================== USER INPUTS ==================

% Time crop
t_start = 0;  % seconds
t_end   = 10000;  % seconds

% Dataset format flag
% 'no_drs'  - 21 columns, no DRS channel
% 'drs'     - 22 columns, DRS but no IMU channels
% 'drs_imu' - 27 columns, DRS and full IMU channels
dataset_format = 'drs_imu';

% Vehicle parameters
VehicleMass     = 259;      % kg
RotorMass_front = 0.5107;   % kg
RotorMass_rear  = 0.3343;   % kg
RotorArea_front = 0.038;    % m^2
RotorArea_rear  = 0.0226;   % m^2
Cd              = 1.46;     % drag coefficient
AeroA           = 1;        % cross sectional area m^2
I               = 0.30754;  % rotational inertia kg*m^2
WheelR          = 0.213;    % meters
AirDen          = 1.14;      % kg/m^3
gear_ratio      = 12.97;    % motor to wheel gear ratio

% Ambient temperature
TambC = 25;  % degC

% Brake geometry
front_piston_count = 6;
front_piston_dia   = 0.0157226;  % meters
front_rotor_dia    = 0.18542;    % meters
rear_piston_count  = 4;
rear_piston_dia    = 0.0141732;  % meters
rear_rotor_dia     = 0.18796;    % meters
front_MC_dia       = 0.01778;    % meters
rear_MC_dia        = 0.01778;    % meters
pedal_ratio        = 4.558;

% Brake temp sim fit parameters
x1 = 2.0006;    % convection coefficient slope
b1 = 40.7137;   % convection coefficient intercept
x2 = 0.001005;  % pad fraction slope
b2 = -0.5385;   % pad fraction intercept

% Mu vs temperature lookup table
mu_temp_table = [100, 200, 300, 400, 500, 600, 700, 800, 900, 950, 1100, 1200];  % degF
mu_table      = [0.45, 0.46, 0.49, 0.53, 0.55, 0.56, 0.565, 0.565, 0.57, 0.57, 0.535, 0.535];

% Lockup torque tables
lockup_speeds      = [5, 15, 25, 35, 45, 55, 65, 75];  % mph
front_lockup_table = [246.8879, 252.0283, 261.8657, 275.4943, 291.7, 309.0181, 325.8681, 340.5926];  % Nm
rear_lockup_table  = [105.8091, 108.0121, 112.2281, 118.069, 125.0143, 132.4363, 139.6577, 145.9683];  % Nm

% Regen parameters
current_limit_actual     = 35;           % A, current setting
current_limits           = [30, 35, 40]; % A, limits to compare
torque_limit_total_MN    = 200;          % %Mn, actual observed maximum from data
torque_limit_total_wheel = torque_limit_total_MN/100*9.8*gear_ratio;  % Nm, total across all 4 wheels
regen_power_limit        = 30000;        % W, motor controller power limit (set from controller or fit from data)
regen_max_curve_a        = 7041;         % coefficient for regen max curve
regen_max_curve_b        = -0.935;       % exponent for regen max curve
regen_efficiency         = 0.66;         % utilization factor to match actual data
min_regen_speed          = 5;            % mph, minimum speed for regen (system limit)
min_energy_threshold     = 50;           % kJ, minimum energy before calculating percent difference

% Aero drag curve coefficients (drag force in N vs speed in mph)
% DRS Open:   F = a*v^2 + b*v + c
% DRS Closed: F = a*v^2 + b*v + c
aero_open_a   =  0.0793024;   aero_open_b   =  1.46483;   aero_open_c   = -37.23993;
aero_closed_a =  0.136247;    aero_closed_b =  2.53449;    aero_closed_c = -67.44975;

% Velocity cleaning threshold
velx_threshold = 34;  % m/s

% Regen index thresholds
regen_torque_threshold = -10;   % Nm wheel side
regen_accel_threshold  = -1.5;  % m/s^2

% Plot colors for current limits
colors = {'r-', 'g-', 'b-'};

%% ================== DATA LOADING ==================
[file, path] = uigetfile('*.txt', 'Select dataset file');
filename = fullfile(path, file);

data = readmatrix(filename, ...
    'FileType', 'text', ...
    'Delimiter', '\t', ...
    'NumHeaderLines', 6, ...
    'TreatAsMissing', {'', 'NaN'});

% Crop by time
t_raw = data(:, 1);
data  = data(t_raw >= t_start & t_raw <= t_end, :);

%% ================== DATA PARSING ==================
% Column order varies by dataset format - see dataset_format flag above
%
% drs_imu (27 columns):
% 1: xtime, 2-5: temp ADC FL/FR/RL/RR, 6: DRS_State, 7: BMS_Current,
% 8: brake_pressure_front, 9: SteerPct, 10-11: AccelPositionA/B,
% 12: BrakesRear_ADC, 13-20: Motor Torque/Velocity FL/FR/RL/RR,
% 21: AccelY, 22: AccelZ, 23-25: AngularRate X/Y/Z, 26: VelBodyX, 27: AccelX
%
% drs (22 columns): same as drs_imu but without columns 21-25
% no_drs (21 columns): same as drs but without column 6 (DRS_State)

t           = data(:, 1);
fl_temp_adc = data(:, 2);
fr_temp_adc = data(:, 3);
rl_temp_adc = data(:, 4);
rr_temp_adc = data(:, 5);

if strcmp(dataset_format, 'drs_imu')
    drs_state         = data(:, 6);
    bms_current       = data(:, 7);
    frontpressure_adc = data(:, 8);
    steer_pct         = data(:, 9);
    accel_posA        = data(:, 10);
    accel_posB        = data(:, 11);
    rearpressure_adc  = data(:, 12);
    fl_Tmotor_Mn      = data(:, 13);
    fl_vmotor         = data(:, 14);
    fr_Tmotor_Mn      = data(:, 15);
    fr_vmotor         = data(:, 16);
    rl_Tmotor_Mn      = data(:, 17);
    rl_vmotor         = data(:, 18);
    rr_Tmotor_Mn      = data(:, 19);
    rr_vmotor         = data(:, 20);
    accely            = data(:, 21);  % m/s^2 lateral
    accelz            = data(:, 22);  % m/s^2 vertical
    angular_rate_x    = data(:, 23);  % rad/s roll rate
    angular_rate_y    = data(:, 24);  % rad/s pitch rate
    angular_rate_z    = data(:, 25);  % rad/s yaw rate
    velx              = data(:, 26);  % m/s
    accelx            = data(:, 27);  % m/s^2
elseif strcmp(dataset_format, 'drs')
    drs_state         = data(:, 6);
    bms_current       = data(:, 7);
    frontpressure_adc = data(:, 8);
    steer_pct         = data(:, 9);
    accel_posA        = data(:, 10);
    accel_posB        = data(:, 11);
    rearpressure_adc  = data(:, 12);
    fl_Tmotor_Mn      = data(:, 13);
    fl_vmotor         = data(:, 14);
    fr_Tmotor_Mn      = data(:, 15);
    fr_vmotor         = data(:, 16);
    rl_Tmotor_Mn      = data(:, 17);
    rl_vmotor         = data(:, 18);
    rr_Tmotor_Mn      = data(:, 19);
    rr_vmotor         = data(:, 20);
    velx              = data(:, 21);
    accelx            = data(:, 22);
    accely            = zeros(size(t));
    accelz            = zeros(size(t));
    angular_rate_x    = zeros(size(t));
    angular_rate_y    = zeros(size(t));
    angular_rate_z    = zeros(size(t));
elseif strcmp(dataset_format, 'no_drs')
    drs_state         = zeros(size(data, 1), 1);  % default DRS closed
    bms_current       = data(:, 6);
    frontpressure_adc = data(:, 7);
    steer_pct         = data(:, 8);
    accel_posA        = data(:, 9);
    accel_posB        = data(:, 10);
    rearpressure_adc  = data(:, 11);
    fl_Tmotor_Mn      = data(:, 12);
    fl_vmotor         = data(:, 13);
    fr_Tmotor_Mn      = data(:, 14);
    fr_vmotor         = data(:, 15);
    rl_Tmotor_Mn      = data(:, 16);
    rl_vmotor         = data(:, 17);
    rr_Tmotor_Mn      = data(:, 18);
    rr_vmotor         = data(:, 19);
    velx              = data(:, 20);
    accelx            = data(:, 21);
    accely            = zeros(size(t));
    accelz            = zeros(size(t));
    angular_rate_x    = zeros(size(t));
    angular_rate_y    = zeros(size(t));
    angular_rate_z    = zeros(size(t));
end

%% ================== DATA CONVERSION ==================
frontpressure = 0.924 * frontpressure_adc - 332.64;  % psi
rearpressure  = 0.924 * rearpressure_adc  - 376.068; % psi
frontpressure = max(frontpressure, 0);
rearpressure  = max(rearpressure,  0);

fl_temp_C = 0.246 * (fl_temp_adc - 406);  % degC
fr_temp_C = 0.246 * (fr_temp_adc - 406);
rl_temp_C = 0.246 * (rl_temp_adc - 406);
rr_temp_C = 0.246 * (rr_temp_adc - 406);
fl_temp_F = fl_temp_C * (9/5) + 32;  % degF
fr_temp_F = fr_temp_C * (9/5) + 32;
rl_temp_F = rl_temp_C * (9/5) + 32;
rr_temp_F = rr_temp_C * (9/5) + 32;

fl_Tmotor = fl_Tmotor_Mn / 100 * 9.8;  % Nm
fr_Tmotor = fr_Tmotor_Mn / 100 * 9.8;
rl_Tmotor = rl_Tmotor_Mn / 100 * 9.8;
rr_Tmotor = rr_Tmotor_Mn / 100 * 9.8;
fl_Twheel = fl_Tmotor * gear_ratio;  % Nm
fr_Twheel = fr_Tmotor * gear_ratio;
rl_Twheel = rl_Tmotor * gear_ratio;
rr_Twheel = rr_Tmotor * gear_ratio;
fl_vwheel = fl_vmotor / gear_ratio;  % rpm
fr_vwheel = fr_vmotor / gear_ratio;
rl_vwheel = rl_vmotor / gear_ratio;
rr_vwheel = rr_vmotor / gear_ratio;

accel_pos   = (accel_posA + accel_posB) / 2;
steer_angle = steer_pct;
pbias       = frontpressure ./ (frontpressure + rearpressure);
pbias(isnan(pbias)) = 0;

velx(abs(velx) > velx_threshold) = 0;
velx(isnan(velx)) = 0;
velx(velx < 0)    = 0;
speed_mph = velx * 2.23694;  % m/s to mph after cleaning

%% ================== COMPUTATIONS ==================

% Derived brake geometry
front_piston_area  = front_piston_count * pi * (front_piston_dia/2)^2;
rear_piston_area   = rear_piston_count  * pi * (rear_piston_dia/2)^2;
front_rotor_radius = front_rotor_dia / 2;
rear_rotor_radius  = rear_rotor_dia  / 2;

% Pedal force
MC_area_front   = pi * (front_MC_dia/2)^2;  % m^2
MC_area_rear    = pi * (rear_MC_dia/2)^2;   % m^2
pedal_force     = (frontpressure*6895.*MC_area_front + rearpressure*6895.*MC_area_rear) / pedal_ratio;  % N
pedal_force_lbf = pedal_force / 4.44822;  % lbf

% Clamp forces
fl_clamp_force = (frontpressure * 6895) .* front_piston_area;
fr_clamp_force = (frontpressure * 6895) .* front_piston_area;
rl_clamp_force = (rearpressure  * 6895) .* rear_piston_area;
rr_clamp_force = (rearpressure  * 6895) .* rear_piston_area;

% Temperature dependent mu
mu_front = interp1(mu_temp_table, mu_table, fr_temp_F, 'linear', 'extrap');
mu_front = max(min(mu_front, max(mu_table)), min(mu_table));
mu_rear  = interp1(mu_temp_table, mu_table, rr_temp_F, 'linear', 'extrap');
mu_rear  = max(min(mu_rear,  max(mu_table)), min(mu_table));

% Brake torques
fl_Tbrake   = -2 * mu_front .* fl_clamp_force .* front_rotor_radius;
fr_Tbrake   = -2 * mu_front .* fr_clamp_force .* front_rotor_radius;
rl_Tbrake   = -2 * mu_rear  .* rl_clamp_force .* rear_rotor_radius;
rr_Tbrake   = -2 * mu_rear  .* rr_clamp_force .* rear_rotor_radius;
ftot_Tbrake = fl_Tbrake + fr_Tbrake;
rtot_Tbrake = rl_Tbrake + rr_Tbrake;
Tbias_brake = ftot_Tbrake ./ (rtot_Tbrake + ftot_Tbrake);
Tbias_brake(isnan(Tbias_brake)) = 0;

% Total torques
fl_Ttotal   = fl_Twheel + fl_Tbrake;
fr_Ttotal   = fr_Twheel + fr_Tbrake;
rl_Ttotal   = rl_Twheel + rl_Tbrake;
rr_Ttotal   = rr_Twheel + rr_Tbrake;
ftot_Ttotal = fl_Ttotal + fr_Ttotal;
rtot_Ttotal = rl_Ttotal + rr_Ttotal;
Tbias_decel = ftot_Ttotal ./ (ftot_Ttotal + rtot_Ttotal);
Tbias_decel(isnan(Tbias_decel)) = 0;

% Lockup check
front_lockup_torque = interp1(lockup_speeds, front_lockup_table, speed_mph, 'linear', 'extrap');
rear_lockup_torque  = interp1(lockup_speeds, rear_lockup_table,  speed_mph, 'linear', 'extrap');
front_lockup_torque = max(min(front_lockup_torque, front_lockup_table(end)), front_lockup_table(1));
rear_lockup_torque  = max(min(rear_lockup_torque,  rear_lockup_table(end)),  rear_lockup_table(1));
fl_lockup_ratio     = fl_Ttotal ./ front_lockup_torque;
fr_lockup_ratio     = fr_Ttotal ./ front_lockup_torque;
rl_lockup_ratio     = rl_Ttotal ./ rear_lockup_torque;
rr_lockup_ratio     = rr_Ttotal ./ rear_lockup_torque;

% Max temperatures
[fl_temp_F_max, fl_temp_F_max_idx] = max(fl_temp_F);
[fr_temp_F_max, fr_temp_F_max_idx] = max(fr_temp_F);
[rl_temp_F_max, rl_temp_F_max_idx] = max(rl_temp_F);
[rr_temp_F_max, rr_temp_F_max_idx] = max(rr_temp_F);
fl_temp_F_max_t = t(fl_temp_F_max_idx);
fr_temp_F_max_t = t(fr_temp_F_max_idx);
rl_temp_F_max_t = t(rl_temp_F_max_idx);
rr_temp_F_max_t = t(rr_temp_F_max_idx);

fprintf('FL Max Temp: %.2f degF at t = %.3f s\n', fl_temp_F_max, fl_temp_F_max_t);
fprintf('FR Max Temp: %.2f degF at t = %.3f s\n', fr_temp_F_max, fr_temp_F_max_t);
fprintf('RL Max Temp: %.2f degF at t = %.3f s\n', rl_temp_F_max, rl_temp_F_max_t);
fprintf('RR Max Temp: %.2f degF at t = %.3f s\n', rr_temp_F_max, rr_temp_F_max_t);

% Time step
dt = diff(t);
dt = [dt(1); dt];

% Angular velocities
fl_omega_wheel  = fl_vwheel * (2*pi/60);
fr_omega_wheel  = fr_vwheel * (2*pi/60);
rl_omega_wheel  = rl_vwheel * (2*pi/60);
rr_omega_wheel  = rr_vwheel * (2*pi/60);
fl_omega_motor  = fl_vmotor * (2*pi/60);
fr_omega_motor  = fr_vmotor * (2*pi/60);
rl_omega_motor  = rl_vmotor * (2*pi/60);
rr_omega_motor  = rr_vmotor * (2*pi/60);

% Aero drag force using DRS-dependent polynomial (speed in mph, force in N)
F_aero_open   = aero_open_a   * speed_mph.^2 + aero_open_b   * speed_mph + aero_open_c;
F_aero_closed = aero_closed_a * speed_mph.^2 + aero_closed_b * speed_mph + aero_closed_c;
F_aero = F_aero_open .* double(drs_state == 1) + F_aero_closed .* double(drs_state == 0);
F_aero = max(F_aero, 0);
P_aero = F_aero .* velx;  % W

% Mechanical braking energy
mech_braking_idx = (frontpressure > 5 | rearpressure > 5);
fl_mech_power    = abs(fl_Tbrake .* fl_omega_wheel) .* double(mech_braking_idx);
fr_mech_power    = abs(fr_Tbrake .* fr_omega_wheel) .* double(mech_braking_idx);
rl_mech_power    = abs(rl_Tbrake .* rl_omega_wheel) .* double(mech_braking_idx);
rr_mech_power    = abs(rr_Tbrake .* rr_omega_wheel) .* double(mech_braking_idx);
total_mech_power = fl_mech_power + fr_mech_power + rl_mech_power + rr_mech_power;
mechanical_energy_J  = sum(total_mech_power .* dt);
mechanical_energy_kJ = mechanical_energy_J / 1000;
mechanical_energy_Wh = mechanical_energy_J / 3600;

% Regenerative braking energy
decelerating_idx = accelx < 0;
fl_regen_power   = min(fl_Tmotor .* fl_omega_motor, 0) .* double(decelerating_idx);
fr_regen_power   = min(fr_Tmotor .* fr_omega_motor, 0) .* double(decelerating_idx);
rl_regen_power   = min(rl_Tmotor .* rl_omega_motor, 0) .* double(decelerating_idx);
rr_regen_power   = min(rr_Tmotor .* rr_omega_motor, 0) .* double(decelerating_idx);
total_regen_power     = fl_regen_power + fr_regen_power + rl_regen_power + rr_regen_power;
total_regen_energy_J  = abs(sum(total_regen_power .* dt));
total_regen_energy_kJ = total_regen_energy_J / 1000;
total_regen_energy_Wh = total_regen_energy_J / 3600;
total_braking_energy_J  = mechanical_energy_J + total_regen_energy_J;
total_braking_energy_kJ = total_braking_energy_J / 1000;
total_braking_energy_Wh = total_braking_energy_J / 3600;
regen_fraction = total_regen_energy_J / total_braking_energy_J * 100;

fprintf('\n--- Braking Energy Summary ---\n');
fprintf('Total Braking Energy:   %.2f kJ  (%.2f Wh)\n', total_braking_energy_kJ, total_braking_energy_Wh);
fprintf('Regenerative Braking:   %.2f kJ  (%.2f Wh)\n', total_regen_energy_kJ,   total_regen_energy_Wh);
fprintf('Mechanical Braking:     %.2f kJ  (%.2f Wh)\n', mechanical_energy_kJ,    mechanical_energy_Wh);
fprintf('Regen Fraction:         %.1f%%\n',              regen_fraction);
fprintf('Number of Braking Events: %d\n',               sum(diff(mech_braking_idx) == 1));

% Regen wheel torques and bias
fl_Twheel_regen    = min(fl_Twheel, 0);
fr_Twheel_regen    = min(fr_Twheel, 0);
rl_Twheel_regen    = min(rl_Twheel, 0);
rr_Twheel_regen    = min(rr_Twheel, 0);
ftot_Twheel_regen  = fl_Twheel_regen + fr_Twheel_regen;
rtot_Twheel_regen  = rl_Twheel_regen + rr_Twheel_regen;
total_Twheel_regen = ftot_Twheel_regen + rtot_Twheel_regen;
Tbias_regen        = ftot_Twheel_regen ./ (ftot_Twheel_regen + rtot_Twheel_regen);
Tbias_regen(isnan(Tbias_regen)) = 0;

% Event indices
regen_idx   = bms_current < 0 & speed_mph > min_regen_speed;
braking_idx = fl_Ttotal < 0 | fr_Ttotal < 0 | rl_Ttotal < 0 | rr_Ttotal < 0;
decel_idx   = accelx < 0;
t_brake     = t(braking_idx);

% Aero drag energy during regen events
aero_energy_during_regen_J  = sum(P_aero .* double(regen_idx) .* dt);
aero_energy_during_regen_kJ = aero_energy_during_regen_J / 1000;

fprintf('\n--- Aero Drag During Regen ---\n');
fprintf('Aero Energy during regen: %.2f kJ\n', aero_energy_during_regen_kJ);

% Cumulative regen energy
cumulative_regen_energy_kJ = cumsum(abs(total_regen_power .* dt)) / 1000;

% Wheel omega from speed for power limit calculation
omega_wheel_from_speed       = speed_mph * 0.44704 / WheelR;  % rad/s
omega_wheel_from_speed       = max(omega_wheel_from_speed, 0.1);  % avoid divide by zero
torque_from_power_limit      = regen_power_limit ./ omega_wheel_from_speed;  % Nm

% Predicted regen at different current limits
cumulative_regen_pred = zeros(length(t), length(current_limits));

for i = 1:length(current_limits)
    scale_factor           = current_limits(i) / current_limit_actual;
    max_regen_torque       = regen_max_curve_a * speed_mph.^(regen_max_curve_b) * scale_factor * regen_efficiency;
    max_regen_torque(speed_mph < 1) = 0;

    % Apply all three limits: current-dependent curve, power limit, torque limit
    % Take the minimum (most restrictive) at each timestep
    effective_max_torque   = min(min(max_regen_torque, torque_from_power_limit), torque_limit_total_wheel);

    predicted_regen_torque = -effective_max_torque .* double(regen_idx);
    predicted_regen_power  = (predicted_regen_torque .* Tbias_regen / 2) .* fl_omega_wheel + ...
                             (predicted_regen_torque .* Tbias_regen / 2) .* fr_omega_wheel + ...
                            (predicted_regen_torque .* (1-Tbias_regen) / 2) .* rl_omega_wheel + ...
                              (predicted_regen_torque .* (1-Tbias_regen) / 2) .* rr_omega_wheel;
    cumulative_regen_pred(:, i) = cumsum(abs(predicted_regen_power .* dt)) / 1000;
end

% Percent difference — suppress until min energy threshold is reached
idx_actual = current_limits == current_limit_actual;
pct_diff   = zeros(size(t));
valid_idx  = cumulative_regen_energy_kJ > min_energy_threshold;
pct_diff(valid_idx) = (cumulative_regen_pred(valid_idx, idx_actual) - cumulative_regen_energy_kJ(valid_idx)) ./ ...
                       cumulative_regen_energy_kJ(valid_idx) * 100;

fprintf('\n--- %dA Regen Model vs Actual ---\n', current_limit_actual);
fprintf('Actual Regen Energy:    %.2f kJ\n', cumulative_regen_energy_kJ(end));
fprintf('Predicted Regen Energy: %.2f kJ\n', cumulative_regen_pred(end, idx_actual));
fprintf('Percent Difference:     %.2f%%\n',  pct_diff(end));
fprintf('Aero fraction of predicted: %.1f%%\n', aero_energy_during_regen_kJ / cumulative_regen_pred(end, idx_actual) * 100);

% Cumulative energy breakdown
cumulative_mech_energy_kJ = cumsum(total_mech_power .* dt) / 1000;
cumulative_aero_energy_kJ = cumsum(P_aero .* double(decel_idx) .* dt) / 1000;
cumulative_total_decel_kJ = cumulative_mech_energy_kJ + cumulative_regen_energy_kJ + cumulative_aero_energy_kJ;
olddrag = 0.0022*(fl_vmotor+fr_vmotor+rl_vmotor+rr_vmotor); %W power loss
newdrag = 0.006*(fl_vmotor+fr_vmotor+rl_vmotor+rr_vmotor); %W power loss
cumulative_brakedrag_old_kJ = cumsum(olddrag.* dt) / 1000;
cumulative_brakedrag_new_kJ = cumsum(newdrag.* dt) / 1000;

fprintf('\n--- Braking Energy Breakdown ---\n');
fprintf('Mechanical:   %.2f kJ  (%.1f%%)\n', cumulative_mech_energy_kJ(end),  cumulative_mech_energy_kJ(end)  / cumulative_total_decel_kJ(end) * 100);
fprintf('Regenerative: %.2f kJ  (%.1f%%)\n', cumulative_regen_energy_kJ(end), cumulative_regen_energy_kJ(end) / cumulative_total_decel_kJ(end) * 100);
fprintf('Aero:         %.2f kJ  (%.1f%%)\n', cumulative_aero_energy_kJ(end),  cumulative_aero_energy_kJ(end)  / cumulative_total_decel_kJ(end) * 100);
fprintf('Total:        %.2f kJ\n',           cumulative_total_decel_kJ(end));

% DRS drag energy analysis
aero_energy_open_kJ       = sum(F_aero_open  .* max(velx, 0) .* double(drs_state == 1) .* dt) / 1000;
aero_energy_closed_kJ     = sum(F_aero_closed .* max(velx, 0) .* double(drs_state == 0) .* dt) / 1000;
aero_energy_all_open_kJ   = sum(F_aero_open   .* max(velx, 0) .* dt) / 1000;
aero_energy_all_closed_kJ = sum(F_aero_closed  .* max(velx, 0) .* dt) / 1000;

fprintf('\n--- DRS Drag Energy Analysis ---\n');
fprintf('Actual aero energy (DRS open):      %.2f kJ\n', aero_energy_open_kJ);
fprintf('Actual aero energy (DRS closed):    %.2f kJ\n', aero_energy_closed_kJ);
fprintf('Total actual aero energy:           %.2f kJ\n', aero_energy_open_kJ + aero_energy_closed_kJ);
fprintf('If always DRS open:                 %.2f kJ\n', aero_energy_all_open_kJ);
fprintf('If always DRS closed:               %.2f kJ\n', aero_energy_all_closed_kJ);
fprintf('DRS energy saving vs always closed: %.2f kJ (%.1f%%)\n', ...
    aero_energy_all_closed_kJ - (aero_energy_open_kJ + aero_energy_closed_kJ), ...
    (aero_energy_all_closed_kJ - (aero_energy_open_kJ + aero_energy_closed_kJ)) / aero_energy_all_closed_kJ * 100);

% % Manual deceleration event analysis
% decel_t_start = 4203.9;  % seconds, start of braking event
% decel_t_end   = 4206.2;  % seconds, end of braking event
% 
% decel_event_idx   = t >= decel_t_start & t <= decel_t_end;
% avg_decel_g       = mean(accelx(decel_event_idx)) / 9.81;
% avg_decel_ms2     = mean(accelx(decel_event_idx));
% peak_decel_g      = min(accelx(decel_event_idx)) / 9.81;
% entry_speed_mph   = speed_mph(find(decel_event_idx, 1, 'first'));
% exit_speed_mph    = speed_mph(find(decel_event_idx, 1, 'last'));
% 
% fprintf('\n--- Manual Decel Event Analysis ---\n');
% fprintf('Time window:       %.2f - %.2f s (%.2f s duration)\n', decel_t_start, decel_t_end, decel_t_end - decel_t_start);
% fprintf('Entry speed:       %.2f mph\n', entry_speed_mph);
% fprintf('Exit speed:        %.2f mph\n', exit_speed_mph);
% fprintf('Avg decel rate:    %.3f g  (%.3f m/s^2)\n', avg_decel_g, avg_decel_ms2);
% fprintf('Peak decel rate:   %.3f g\n', peak_decel_g);

% Brake temp sim
[sim_temp_F_front, sim_temp_F_rear, sim_err_front, sim_err_rear] = ...
    brake_temp_sim(t, velx, frontpressure, rearpressure, fr_temp_F, rr_temp_F, ...
    VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
    Tbias_brake, I, WheelR, TambC, x1, b1, x2, b2, ...
    total_regen_power, F_aero, fl_omega_wheel, fr_omega_wheel, rl_omega_wheel, rr_omega_wheel);

%% ================== PLOTS ==================
corners       = {'FL', 'FR', 'RL', 'RR'};
wheel_torques = {fl_Twheel, fr_Twheel, rl_Twheel, rr_Twheel};
brake_torques = {fl_Tbrake, fr_Tbrake, rl_Tbrake, rr_Tbrake};
total_torques = {fl_Ttotal, fr_Ttotal, rl_Ttotal, rr_Ttotal};
lockup_ratios = {fl_lockup_ratio, fr_lockup_ratio, rl_lockup_ratio, rr_lockup_ratio};
speed_plot    = linspace(1, 80, 500);  % mph for theoretical curves

input('Press Enter to generate plots...');

%% Plot 1: Brake Temperatures over Time
figure('Name', 'Brake Temps over Time');
plot(t, fr_temp_F, 'b-', 'DisplayName', 'FR'); hold on;
plot(t, rr_temp_F, 'm-', 'DisplayName', 'RR');
xlabel('Time (s)'); ylabel('Temperature (°F)');
title('Brake Temperatures over Time');
legend('Location', 'best'); grid on;

%% Plot 2a: Torques per Corner
figure('Name', 'Torques Per Corner');
for i = 1:4
    subplot(4, 1, i);
    plot(t, wheel_torques{i}, 'b-', 'DisplayName', 'Motor Torque'); hold on;
    plot(t, brake_torques{i}, 'r-', 'DisplayName', 'Brake Torque');
    plot(t, total_torques{i}, '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Total Torque');
    ylim([min([fl_Ttotal; fr_Ttotal; rl_Ttotal; rr_Ttotal]), 0]);
    xlabel('Time (s)'); ylabel('Torque (Nm)');
    title([corners{i} ' Torque over Time']);
    legend('Location', 'best'); grid on;
end
sgtitle('Wheel Torques over Time');

%% Plot 2b: Lockup Ratios per Corner
figure('Name', 'Lockup Ratios');
for i = 1:4
    subplot(4, 1, i);
    plot(t, lockup_ratios{i}, 'k-'); hold on;
    yline(-1.0, 'r--', 'Lockup Threshold');
    xlabel('Time (s)'); ylabel('Lockup Ratio');
    ylim([min([fl_lockup_ratio; fr_lockup_ratio; rl_lockup_ratio; rr_lockup_ratio; -1.2]), 0]);
    title([corners{i} ' Lockup Ratio over Time']); grid on;
end
sgtitle('Wheel Lockup Ratios over Time');

%% Plot 3: Brake Torque Bias over Time
figure('Name', 'Brake Torque Bias');
scatter(t_brake, Tbias_brake(braking_idx), 50, 'b.');
xlabel('Time (s)'); ylabel('Brake Torque Bias');
title('Brake Torque Bias during Braking Events (Front/Total)'); grid on;

%% Plot 4: Total Decel Torque Bias over Time
figure('Name', 'Total Decel Torque Bias');
scatter(t_brake, Tbias_decel(braking_idx), 50, 'k.');
xlabel('Time (s)'); ylabel('Total Decel Torque Bias (Front/Total)');
title('Total Decel Torque Bias during Braking Events'); grid on;

%% Plot 5: Regen Torque Bias vs Speed
figure('Name', 'Regen Torque Bias vs Speed');
scatter(speed_mph(regen_idx), Tbias_regen(regen_idx), 50, 'k.');
xlabel('Vehicle Speed (mph)'); ylabel('Regen Torque Bias (Front/Total)');
title('Regen Torque Bias vs Vehicle Speed'); grid on;

%% Plot 6: Total Decel Torque Bias vs Speed
figure('Name', 'Total Decel Torque Bias vs Speed');

% Only plot where total torque is large enough to give meaningful bias
min_torque_threshold = 20;  % Nm, adjust as needed
min_torque_threshold = 150;  % Nm, increase to filter noise-level events
valid_bias_idx = braking_idx & ...
                 abs(ftot_Ttotal + rtot_Ttotal) > min_torque_threshold & ...
                 speed_mph > 2 & ...
                 accelx < -1;

scatter(speed_mph(valid_bias_idx), Tbias_decel(valid_bias_idx), 50, 'k.');
ylim([0 1]);
xlabel('Vehicle Speed (mph)'); ylabel('Total Decel Torque Bias (Front/Total)');
title('Total Decel Torque Bias vs Vehicle Speed during Braking Events');
grid on;

%% Plot 7: Brake Torque Bias vs Front Brake Pressure
figure('Name', 'Brake Torque Bias vs Front Brake Pressure');
fit_idx = frontpressure > 10 & ~isnan(Tbias_brake);
x_data  = frontpressure(fit_idx);
y_data  = Tbias_brake(fit_idx);
x_fit   = linspace(min(x_data), max(x_data), 100);
p_log   = polyfit(log(x_data), y_data, 1);
y_log   = p_log(1) * log(x_fit) + p_log(2);
scatter(x_data, y_data, 20, 'b.', 'DisplayName', 'Torque Bias'); hold on;
% plot(x_fit, y_log, 'r-', 'LineWidth', 2, 'DisplayName', ...
%     sprintf('Log Fit: y = %.4f*ln(x) + %.4f', p_log(1), p_log(2)));
xlabel('Front Brake Pressure (psi)'); ylabel('Brake Torque Bias (Front/Total)');
title('Brake Torque Bias vs Front Brake Pressure');
legend('Location', 'best'); grid on;

%% Plot 8: Regen Total Wheel Torque vs Speed (colored by BMS current)
figure('Name', 'Regen Total Wheel Torque vs Speed');
scatter(speed_mph(regen_idx), abs(total_Twheel_regen(regen_idx)), 20, bms_current(regen_idx), '.');
colorbar; colormap('jet');
clabel = colorbar; clabel.Label.String = 'Current (A)';
xlabel('Vehicle Speed (mph)'); ylabel('Regen Wheel Torque (Nm)');
title('Regen Total Wheel Torque vs Speed'); grid on;

%% Plot 9: Total Regen Wheel Torque over Time
figure('Name', 'Total Regen Wheel Torque over Time');
plot(t, abs(total_Twheel_regen), 'k-');
xlabel('Time (s)'); ylabel('Total Regen Wheel Torque (Nm)');
title('Total Regen Wheel Torque over Time'); grid on;

%% Plot 10: Simulated vs Measured Brake Temperatures
figure('Name', 'Simulated vs Measured Brake Temps');
subplot(2, 1, 1);
plot(t, fr_temp_F, 'b-', 'DisplayName', 'Measured'); hold on;
plot(t, sim_temp_F_front, 'r-', 'DisplayName', 'Simulated Front');
xlabel('Time (s)'); ylabel('Temperature (°F)');
title('Front Brake Temperature: Simulated vs Measured');
legend('Location', 'best'); grid on;
subplot(2, 1, 2);
plot(t, rr_temp_F, 'b-', 'DisplayName', 'Measured'); hold on;
plot(t, sim_temp_F_rear, 'r-', 'DisplayName', 'Simulated Rear');
xlabel('Time (s)'); ylabel('Temperature (°F)');
title('Rear Brake Temperature: Simulated vs Measured');
legend('Location', 'best'); grid on;
sgtitle('Brake Temperature Simulation vs Measured');

%% Plot 11: Driver Inputs over Time
figure('Name', 'Driver Inputs');
subplot(4, 1, 1);
plot(t, frontpressure, 'b-', 'DisplayName', 'Front'); hold on;
plot(t, rearpressure,  'r-', 'DisplayName', 'Rear');
xlabel('Time (s)'); ylabel('Brake Pressure (psi)');
title('Brake Pressure'); legend('Location', 'best'); grid on;
subplot(4, 1, 2);
plot(t(mech_braking_idx), pedal_force_lbf(mech_braking_idx), 'b.');
xlabel('Time (s)'); ylabel('Pedal Force (lbf)');
title('Brake Pedal Force'); grid on;
subplot(4, 1, 3);
plot(t, accel_pos, 'g-');
xlabel('Time (s)'); ylabel('Throttle Position (%)');
title('Throttle Position'); grid on;
subplot(4, 1, 4);
plot(t, steer_angle, 'k-');
xlabel('Time (s)'); ylabel('Steering Angle (%)');
title('Steering Angle'); grid on;
sgtitle('Driver Inputs over Time');

%% Plot 12: BMS Current over Time
figure('Name', 'BMS Current');
plot(t, bms_current, 'm-');
% ylim([min(bms_current), 0]);
xlabel('Time (s)'); ylabel('Current (A)');
title('BMS Current over Time'); grid on;

%% Plot 13: Cumulative Regen Energy over Time
figure('Name', 'Cumulative Regen Energy');
plot(t, cumulative_regen_energy_kJ, 'b-');
xlabel('Time (s)'); ylabel('Cumulative Regen Energy (kJ)');
title('Cumulative Regenerative Braking Energy over Time'); grid on;

%% Plot 14: Predicted Cumulative Regen Energy vs Current Limit
figure('Name', 'Cumulative Regen Energy vs Current Limit Estimation');
for i = 1:length(current_limits)
    plot(t, cumulative_regen_pred(:, i), colors{i}, ...
        'DisplayName', sprintf('%dA limit', current_limits(i))); hold on;
end
xlabel('Time (s)'); ylabel('Cumulative Regen Energy (kJ)');
title('Predicted Cumulative Regen Energy vs Current Limit');
legend('Location', 'best'); grid on;

%% Plot 15: Actual vs Predicted Cumulative Regen Energy
figure('Name', sprintf('Actual vs Predicted Regen Energy %dA', current_limit_actual));
subplot(2, 1, 1);
plot(t, cumulative_regen_energy_kJ,           'b-', 'DisplayName', 'Actual'); hold on;
plot(t, cumulative_regen_pred(:, idx_actual),  'r-', 'DisplayName', ...
    sprintf('Predicted (%dA)', current_limit_actual));
xlabel('Time (s)'); ylabel('Cumulative Regen Energy (kJ)');
title(sprintf('Actual vs Predicted Cumulative Regen Energy (%dA)', current_limit_actual));
legend('Location', 'best'); grid on;
subplot(2, 1, 2);
plot(t, pct_diff, 'k-');
yline(0, 'r--');
xlabel('Time (s)'); ylabel('Percent Difference (%)');
title('Predicted vs Actual Percent Difference'); grid on;
sgtitle(sprintf('%dA Regen Model Validation', current_limit_actual));

%% Plot 16: Theoretical Regen Torque vs Speed at Different Current Limits
figure('Name', 'Theoretical Regen Torque vs Speed');
omega_plot         = speed_plot * 0.44704 / WheelR;  % rad/s from mph
omega_plot         = max(omega_plot, 0.1);
torque_power_plot  = regen_power_limit ./ omega_plot;  % power limit curve Nm

for i = 1:length(current_limits)
    scale_factor    = current_limits(i) / current_limit_actual;
    max_regen_speed = regen_max_curve_a * speed_plot.^(regen_max_curve_b) * scale_factor;
    effective_limit = min(min(max_regen_speed, torque_power_plot), torque_limit_total_wheel);
    plot(speed_plot, effective_limit, colors{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('%dA limit', current_limits(i))); hold on;
end
plot(speed_plot, torque_power_plot, 'c--', 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('Power Limit (%.0f kW)', regen_power_limit/1000));
scatter(speed_mph(regen_idx), abs(total_Twheel_regen(regen_idx)), 8, ...
    [0.5, 0.5, 0.5], '.', 'DisplayName', 'Measured Regen Torque');
yline(torque_limit_total_wheel, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Torque Limit');
% xline(40, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Speed Limit');
xlabel('Vehicle Speed (mph)'); ylabel('Max Regen Torque (Nm)');
title('Theoretical Max Regen Torque vs Speed at Different Current Limits');
legend('Location', 'best'); grid on;
ylim([0, torque_limit_total_wheel * 1.2]);

%% Plot 17: Cumulative Aero Drag Energy by DRS State
P_aero_open_actual        = F_aero_open  .* max(velx, 0) .* double(drs_state == 1);
P_aero_closed_actual      = F_aero_closed .* max(velx, 0) .* double(drs_state == 0);
cumulative_aero_open_kJ   = cumsum(P_aero_open_actual   .* dt) / 1000;
cumulative_aero_closed_kJ = cumsum(P_aero_closed_actual .* dt) / 1000;
cumulative_aero_total_kJ  = cumulative_aero_open_kJ + cumulative_aero_closed_kJ;
cumulative_aero_all_open_kJ   = cumsum(F_aero_open  .* max(velx, 0) .* dt) / 1000;
cumulative_aero_all_closed_kJ = cumsum(F_aero_closed .* max(velx, 0) .* dt) / 1000;

figure('Name', 'Cumulative Aero Drag Energy by DRS State');
plot(t, cumulative_aero_total_kJ,      'k-',  'DisplayName', 'Actual (mixed DRS)'); hold on;
plot(t, cumulative_aero_all_open_kJ,   'r--', 'DisplayName', 'If always open');
plot(t, cumulative_aero_all_closed_kJ, 'b--', 'DisplayName', 'If always closed');
xlabel('Time (s)'); ylabel('Cumulative Aero Energy (kJ)');
title('Cumulative Aero Drag Energy by DRS State');
legend('Location', 'best'); grid on;

%% Plot 18: Cumulative Braking Energy Breakdown
figure('Name', 'Cumulative Braking Energy Breakdown');
plot(t, cumulative_total_decel_kJ,  'k-',  'LineWidth', 1.5, 'DisplayName', 'Total Decel'); hold on;
plot(t, cumulative_mech_energy_kJ,  'r-',  'DisplayName', 'Mechanical Braking');
plot(t, cumulative_regen_energy_kJ, 'b-',  'DisplayName', 'Regenerative Braking');
plot(t, cumulative_aero_energy_kJ,  'g-',  'DisplayName', 'Aero Drag');
xlabel('Time (s)'); ylabel('Cumulative Energy (kJ)');
title('Cumulative Braking Energy Breakdown over Time');
legend('Location', 'best'); grid on;

%% Plot 19: Steer Angle vs Yaw Rate
figure('Name', 'Steer Angle vs Yaw Rate');
scatter(steer_angle, angular_rate_z, 10, 'k.');
xlabel('Steering Angle (%)'); ylabel('Yaw Rate (rad/s)');
title('Steering Angle vs Yaw Rate'); grid on;

%% Plot 20: Throttle vs Steer Angle vs Yaw Rate Heatmap
figure('Name', 'Throttle vs Steer Angle vs Yaw Rate');

% Create regular grid
throttle_grid    = linspace(min(accel_pos),   max(accel_pos),   50);
steer_grid       = linspace(min(steer_angle), max(steer_angle), 50);
[T_grid, S_grid] = meshgrid(throttle_grid, steer_grid);

% Interpolate scattered data onto grid
Z_grid = griddata(accel_pos, steer_angle, angular_rate_z, T_grid, S_grid, 'linear');

% Smooth the surface
Z_grid = imgaussfilt(Z_grid, 1.5);

% Plot heatmap with contour lines and data point overlay
contourf(T_grid, S_grid, Z_grid, 20, 'LineColor', 'none');
hold on;
contour(T_grid, S_grid, Z_grid, 10, 'LineColor', 'k', 'LineWidth', 0.5);
scatter(accel_pos, steer_angle, 3, 'w.', 'MarkerEdgeAlpha', 0.2, 'DisplayName', 'Data Points');
colormap('jet');
clabel_h = colorbar;
clabel_h.Label.String = 'Yaw Rate (rad/s)';
xlabel('Throttle Position (%)');
ylabel('Steering Angle (%)');
title('Throttle vs Steering Angle vs Yaw Rate');

%% Plot 21: Vehicle Speed vs Steer Angle vs Yaw Rate Heatmap
figure('Name', 'Vehicle Speed vs Steer Angle vs Yaw Rate');

% Create regular grid
speed_grid       = linspace(min(speed_mph),   max(speed_mph),   50);
steer_grid       = linspace(min(steer_angle), max(steer_angle), 50);
[V_grid, S_grid] = meshgrid(speed_grid, steer_grid);

% Interpolate scattered data onto grid
Z_grid = griddata(speed_mph, steer_angle, angular_rate_z, V_grid, S_grid, 'linear');

% Smooth the surface
Z_grid = imgaussfilt(Z_grid, 1.5);

% Plot heatmap with contour lines and data point overlay
contourf(V_grid, S_grid, Z_grid, 20, 'LineColor', 'none');
hold on;
contour(V_grid, S_grid, Z_grid, 10, 'LineColor', 'k', 'LineWidth', 0.5);
scatter(speed_mph, steer_angle, 3, 'w.', 'MarkerEdgeAlpha', 0.2, 'DisplayName', 'Data Points');
colormap('jet');
clabel_h = colorbar;
clabel_h.Label.String = 'Yaw Rate (rad/s)';
xlabel('Vehicle Speed (mph)');
ylabel('Steering Angle (%)');
title('Vehicle Speed vs Steering Angle vs Yaw Rate');

%% Plot 22: Throttle Position vs Longitudinal Acceleration
figure('Name', 'Throttle Position vs Acceleration');
scatter(accel_pos, accelx, 5, 'k.');
xlabel('Throttle Position (%)'); ylabel('Longitudinal Acceleration (m/s^2)');
title('Throttle Position vs Longitudinal Acceleration'); grid on;
%% Plot: Vehicle Velocity and Acceleration
figure('Name', 'Vehicle Velocity and Acceleration');
subplot(3, 1, 1);
plot(t, speed_mph, 'b-');
xlabel('Time (s)'); ylabel('Vehicle Speed (mph)');
title('Vehicle Speed over Time'); grid on;

subplot(3, 1, 2);
plot(t, accelx, 'r-');
yline(0, 'k--');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
title('Longitudinal Acceleration over Time'); grid on;

subplot(3, 1, 3);
plot(t, accelx / 9.81, 'r-');
yline(0, 'k--');
ylabel('Acceleration (g)');
xlabel('Time (s)');
title('Longitudinal Acceleration over Time'); grid on;

sgtitle('Vehicle Velocity and Acceleration');




% Total decel force minus aero contribution
brake_force_total = (VehicleMass .* abs(accelx)) - F_aero;  % N
brake_force_total = max(brake_force_total, 0);  % clamp to zero

% Back-calculate mu from front brake force
mu_measured = (brake_force_total .* WheelR .* Tbias_brake) ./ ...
    (2 .* fl_clamp_force .* front_rotor_radius);
%% Plot: Measured Pad Mu vs Temperature
mu_plot_idx = mech_braking_idx & ...
              frontpressure > 100 & ...       % higher pressure threshold
              rearpressure  > 70  & ...       % require some rear pressure too
              brake_force_total > 500 & ...  % require meaningful brake force
              ~isnan(mu_measured) & ...
              mu_measured > 0.3 & ...        % physically reasonable range
              mu_measured < 0.8 & ...        % physically reasonable range
              speed_mph > 5;

figure('Name', 'Measured Pad Mu vs Temperature');
scatter(fr_temp_F(mu_plot_idx), mu_measured(mu_plot_idx), 10, 'b.');
xlabel('Brake Temperature (°F)');
ylabel('Pad Mu');
title('Measured Pad Mu vs Brake Temperature');
ylim([0 1]);
grid on;


%% Plot: Measured Pad Mu vs Temperature and Front Brake Pressure
figure('Name', 'Measured Pad Mu vs Temperature and Pressure');

% Create regular grid
temp_grid     = linspace(min(fr_temp_F(mu_plot_idx)),    max(fr_temp_F(mu_plot_idx)),    40);
pressure_grid = linspace(min(frontpressure(mu_plot_idx)), max(frontpressure(mu_plot_idx)), 40);
[T_grid, P_grid] = meshgrid(temp_grid, pressure_grid);

% Interpolate scattered data onto grid
Z_grid = griddata(fr_temp_F(mu_plot_idx), frontpressure(mu_plot_idx), mu_measured(mu_plot_idx), ...
    T_grid, P_grid, 'linear');

% Smooth the surface
Z_grid = imgaussfilt(Z_grid, 2.5);

% Plot
contourf(T_grid, P_grid, Z_grid, 20, 'LineColor', 'none'); hold on;
contour(T_grid,  P_grid, Z_grid, 10, 'LineColor', 'k', 'LineWidth', 0.5);
scatter(fr_temp_F(mu_plot_idx), frontpressure(mu_plot_idx), 15, 'k.');
colormap('jet');
caxis([0.3 0.65]);
clabel_h = colorbar;
clabel_h.Label.String = 'Pad Mu';
xlabel('Brake Temperature (°F)');
ylabel('Front Brake Pressure (psi)');
title('Measured Pad Mu vs Temperature and Front Brake Pressure');
grid on;

%% Plot 13: Cumulative Brake Drag Energy over Time with Percent Difference
figure('Name', 'Cumulative Brake Drag Energy');

plot(t, cumulative_brakedrag_old_kJ, 'b-', t, cumulative_brakedrag_new_kJ, 'k-');
xlabel('Time (s)'); ylabel('Cumulative Brake Drag Energy (kJ)');
title('Cumulative Brake Drag Energy over Time'); grid on;
legend('F34 (new)', 'F33 (old)', 'Location', 'best');