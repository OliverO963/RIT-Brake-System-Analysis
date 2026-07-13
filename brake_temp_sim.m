function [RotorTempArrayF_front, RotorTempArrayF_rear, AvgPercErr_front, AvgPercErr_rear] = brake_temp_sim(t, velx, BrakePress_front, BrakePress_rear, brakeTempArray_front, brakeTempArray_rear, ...
    VehicleMass, RotorMass_front, RotorMass_rear, RotorArea_front, RotorArea_rear, ...
    Tbias, I, WheelR, TambC, x1, b1, x2, b2, regen_power, F_aero, fl_omega_wheel, fr_omega_wheel, rl_omega_wheel, rr_omega_wheel)


TambK = TambC + 273;

%% Calculate Displacement Array
d(1) = 0;
for n = 2:length(t)
    timestep = t(n) - t(n-1);
    avgspd   = (velx(n-1) + velx(n)) / 2;
    d(n)     = avgspd * timestep;
end

%% Calculate Aero Braking Energy using DRS-dependent drag force
Edrag = zeros(size(t));
for k = 2:length(t)
    Edrag(k) = F_aero(k) * d(k);  % N * m = J
end

%% Run Simulation for Front and Rear
[RotorTempArrayF_front, AvgPercErr_front] = run_sim(t, velx, BrakePress_front, brakeTempArray_front, ...
    VehicleMass, RotorMass_front, RotorArea_front, Tbias,      I, WheelR, Edrag, TambK, x1, b1, x2, b2, regen_power, fl_omega_wheel, fr_omega_wheel);

[RotorTempArrayF_rear, AvgPercErr_rear]   = run_sim(t, velx, BrakePress_rear,  brakeTempArray_rear,  ...
    VehicleMass, RotorMass_rear,  RotorArea_rear,  (1-Tbias), I, WheelR, Edrag, TambK, x1, b1, x2, b2, regen_power, rl_omega_wheel, rr_omega_wheel);

fprintf('Average Brake Temp Sim Error (Front): %.2f%%\n', AvgPercErr_front);
fprintf('Average Brake Temp Sim Error (Rear):  %.2f%%\n', AvgPercErr_rear);

end

%% Internal Helper Function
function [RotorTempArrayF, AvgPercErr] = run_sim(t, velx, BrakePress, brakeTempArray, ...
    VehicleMass, RotorMass, RotorArea, Tbias, I, WheelR, Edrag, TambK, x1, b1, x2, b2, total_regen_power, omega_wheel_L, omega_wheel_R)

Energy1         = zeros(size(t));
Energy2         = zeros(size(t));
Energy          = zeros(size(t));
CorrectedEnergy = zeros(size(t));
Power           = zeros(size(t));
RotorTempArrayK(1) = (brakeTempArray(1) - 32) * (5/9) + 273.15; % degF to K
min_pressure = 5;  % psi, threshold to consider brakes applied

for i = 2:length(t)
    prevSpeed = velx(i-1);
    newSpeed  = velx(i);
    % omegaP    = prevSpeed / WheelR;
    % omegaN    = newSpeed  / WheelR;
    DS        = newSpeed - prevSpeed;
    prevTemp  = RotorTempArrayK(i-1);
    tbrake    = t(i) - t(i-1);

    h_w      = x1 * velx(i) + b1;
    Rrotor   = 1 / (h_w * RotorArea);
    % PadFrac  = prevTemp * x2 + b2;
    PadFrac  = max(min(prevTemp * x2 + b2, 1), 0);
    SpecHeat = (0.0005 * prevTemp + 0.2813) * 1000;

    if DS < -2.5
        DS = 0;  % skip - exceeds 2.5g deceleration limit
    end

    if DS < 0 && BrakePress(i) > min_pressure
        Energy1(i)         = 0.5 * VehicleMass * (prevSpeed^2 - newSpeed^2);
        omegaP             = (omega_wheel_L(i-1) + omega_wheel_R(i-1)) / 2;  % average of two corners
        omegaN             = (omega_wheel_L(i)   + omega_wheel_R(i))   / 2;
        Energy2(i)         = 4 * (0.5 * I * (omegaP^2 - omegaN^2));
        Energy(i)          = Energy1(i) + Energy2(i);
        regen_energy       = abs(total_regen_power(i))*tbrake;
        friction_energy    = max(Energy(i)-regen_energy, 0); % filter for noise
        AeroFrac = min(Edrag(i) / max(Energy(i), 1), 1); % filter for noise
        BrakeFrac          = 0.82;
        CorrectedEnergy(i) = friction_energy * 0.5 * Tbias(i) * (1 - AeroFrac) * (1 - PadFrac) * BrakeFrac * 0.8;
        deltaTK            = CorrectedEnergy(i) / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = deltaTK + prevTemp;
        Power(i)           = CorrectedEnergy(i) / tbrake;
        qout               = (RotorTempArrayK(i) - TambK) / Rrotor;
        Eout               = qout * tbrake;
        deltaTKout         = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = RotorTempArrayK(i) - deltaTKout;

    else
        qout               = (prevTemp - TambK) / Rrotor;
        Eout               = qout * tbrake;
        deltaTKout         = Eout / (RotorMass * SpecHeat);
        RotorTempArrayK(i) = prevTemp - deltaTKout;
        Power(i)           = 0;
    end
end

RotorTempArrayF      = ((RotorTempArrayK - 273.15) * (9/5)) + 32;


% Percent error vs measured
PercErr = zeros(length(RotorTempArrayF), 1);
for n = 1:length(RotorTempArrayF)
    PercErr(n) = (abs(RotorTempArrayF(n) - brakeTempArray(n)) / brakeTempArray(n)) * 100;
    if PercErr(n) > 100
        PercErr(n) = PercErr(n-1);
    end
end
AvgPercErr = mean(PercErr);

end