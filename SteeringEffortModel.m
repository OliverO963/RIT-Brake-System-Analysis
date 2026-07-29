%% SteeringEffortModel.m
% Computes the steering-shaft torque (and driver hand force) needed to
% reach a given tire steer angle, from a PAC2002/MF5.2 tire (.tir) file
% and a steering system's geometry.
%
% METHOD
%   1. Tire self-aligning moment Mz(delta) is evaluated with the pure-slip
%      Magic Formula (PAC2002), using slip angle = commanded tire steer
%      angle. This is the standard quasi-static / low-speed approximation
%      (how Mz curves are measured on a flat-trac rig) -- NOT a moving-
%      vehicle cornering sim (no sideslip/yaw-rate/speed effects).
%   2. Steering geometry maps tire angle -> steering-wheel angle via an
%      idealized crank: rack travel x = armLength*sin(delta), then
%      pinion/bevel-gear/column ratios on top of that. The LOCAL ratio
%      d(wheel angle)/d(tire angle) is computed numerically (it's not
%      constant, because the crank kinematics are nonlinear).
%   3. Shaft torque follows from energy conservation (ideal, lossless
%      mechanism): T_wheel * dTheta_wheel = numTires * Mz * dDelta, i.e.
%      T_wheel = numTires * Mz / ratio. Real friction losses would only
%      increase effort above this, so treat results as a lower bound.
%   Not included: Ackermann angle difference / load transfer between the
%   two tires, overturning moment (Mx), scrub radius / kingpin offset.
%
% >>> EDIT THE "USER INPUT" SECTION BELOW <<<

clear; clc; close all;

%% ------------------------- USER INPUT -----------------------------
tirFile = "G:\Shared drives\RIT Formula SAE\Knowledge Center\Vehicle Dynamics\In-House VD Tools\Car Goals Models\Tire Models\AgileTireR20_AdjustedFX.tir";  % path to your PAC2002/MF5.2 .tir file

% VD-driven parameters
geom.armLength           = 0.0853;   % [m]     steering arm length (kingpin to tie-rod ball joint)
geom.rackTravel          = 0.0407; % [m] steering rack travel (0 degrees to full lock in one direction)
geom.maxTireRotationDeg  = (29.963 + 29.484)/2;      % [deg]   max tire steer angle from center, one direction (average of inner and outer to account for anti-Ackermann

% BDC-driven parameters
geom.maxWheelRotationDeg = 100;     % [deg]   max steering-wheel rotation from center, one direction
geom.rackRatio           = geom.rackTravel/(geom.maxWheelRotationDeg/360);    % [m/rev] rack travel per ONE FULL revolution of the pinion
geom.bevelGearRatio      = 1.0;    % [-]     column rotation / pinion rotation (1 = no bevel box)
geom.wheelRadius         = 0.115;    % [m]     steering wheel rim radius

geom.Fz        = [];  % [N]   vertical load per tire; [] -> use tire's nominal load (LFZ0*FNOMIN)
geom.camberDeg = -1.5;   % [deg] static camber angle
geom.numTires  = 2;   % number of steered tires reacting torque through this shaft

nPoints = 200;        % resolution of the sweep from 0 to maxTireRotationDeg
%% ---------------------------------------------------------------------

P = readTIR(tirFile);
if isempty(geom.Fz)
    geom.Fz = P.LFZ0 * P.FNOMIN;
end

tireAngleSweepDeg = linspace(0.05, geom.maxTireRotationDeg, nPoints);
result = steeringShaftTorque(tireAngleSweepDeg, P, geom);

fprintf('--- Geometry consistency check ---\n');
fprintf('At max tire angle (%.1f deg), this geometry gives a wheel angle of %.1f deg ', ...
        geom.maxTireRotationDeg, result.check.predictedWheelAngleAtMaxTire_deg);
fprintf('(you specified %.1f deg, %.1f%% difference).\n', ...
        result.check.specifiedMaxWheelRotationDeg, result.check.percentDifference);
fprintf('If this is off, adjust armLength / rackRatio / bevelGearRatio.\n\n');

[peakTorque, idxPeak] = max(result.steeringTorque_Nm);
fprintf('Peak steering-shaft torque: %.2f N*m at tire angle %.1f deg (wheel angle %.1f deg), ', ...
        peakTorque, result.tireAngleDeg(idxPeak), result.steeringWheelAngleDeg(idxPeak));
fprintf('hand force %.1f N.\n\n', result.handForce_N(idxPeak));

%% ------------------------------- Plots --------------------------------
figure('Name', 'Steering Effort', 'Color', 'w', 'Position', [100 100 1000 700]);

subplot(2,2,1);
plot(result.tireAngleDeg, result.MzPerTire_Nm, 'LineWidth', 1.5); grid on;
xlabel('Tire steer angle [deg]'); ylabel('|M_z| per tire [N\cdotm]');
title('Tire self-aligning moment (Pacejka Mz_0)');

subplot(2,2,2);
plot(result.tireAngleDeg, result.steeringWheelAngleDeg, 'LineWidth', 1.5); grid on;
xlabel('Tire steer angle [deg]'); ylabel('Steering wheel angle [deg]');
title('Steering system kinematics');

subplot(2,2,3);
plot(result.steeringWheelAngleDeg, result.steeringTorque_Nm, 'LineWidth', 1.5); grid on;
xlabel('Steering wheel angle [deg]'); ylabel('Steering shaft torque [N\cdotm]');
title('Torque driver must react, vs. wheel angle');

subplot(2,2,4);
plot(result.steeringWheelAngleDeg, result.handForce_N, 'LineWidth', 1.5); grid on;
xlabel('Steering wheel angle [deg]'); ylabel('Equivalent rim hand force [N]');
title(sprintf('Hand force (rim radius = %.0f mm)', geom.wheelRadius*1000));

sgtitle('Steering Effort Model');


%% ======================================================================
%  LOCAL FUNCTIONS (script-local functions, MATLAB R2016b+)
%  ======================================================================

function out = steeringShaftTorque(tireAngleDeg, tirInput, geom)
%STEERINGSHAFTTORQUE Steering-shaft torque to hold/achieve a tire angle.
%   OUT = STEERINGSHAFTTORQUE(TIREANGLEDEG_deg, TIRINPUT, GEOM) where
%   TIRINPUT is a .tir path or a struct from READTIR, and GEOM has fields
%   armLength, rackRatio, bevelGearRatio, wheelRadius, maxWheelRotationDeg,
%   maxTireRotationDeg, and optional Fz / camberDeg / numTires.
%   Returns OUT.tireAngleDeg, .steeringWheelAngleDeg, .steeringRatio,
%   .MzPerTire_Nm, .steeringTorque_Nm, .handForce_N, and OUT.check (a
%   sanity check comparing the geometry's implied wheel angle at
%   maxTireRotationDeg against the specified maxWheelRotationDeg).

    if ischar(tirInput) || isstring(tirInput)
        P = readTIR(char(tirInput));
    else
        P = tirInput;
    end
    if ~isfield(geom,'Fz') || isempty(geom.Fz),               geom.Fz = P.LFZ0*P.FNOMIN; end
    if ~isfield(geom,'camberDeg') || isempty(geom.camberDeg), geom.camberDeg = 0;        end
    if ~isfield(geom,'numTires') || isempty(geom.numTires),   geom.numTires = 2;         end

    tireAngleDeg  = tireAngleDeg(:);
    deltaTire_rad = deg2rad(tireAngleDeg);
    gammaRad      = deg2rad(geom.camberDeg) * ones(size(deltaTire_rad));

    % Tire self-aligning moment at the commanded angle (magnitude matters
    % for effort, not sign/direction).
    [~, Mz0] = pacejkaPureLateral(deltaTire_rad, geom.Fz, gammaRad, P);
    MzPerTire = abs(Mz0);

    % Steering system kinematics: tire angle -> wheel angle + local ratio
    [thetaSW_rad, ratio] = steeringGeometryAngle(deltaTire_rad, geom);
    steeringWheelAngleDeg = rad2deg(thetaSW_rad);
    ratioMag = max(abs(ratio), 1e-6); % guard divide-by-zero near delta=0

    % Virtual work: T_wheel * dTheta_wheel = numTires * Mz * dDelta
    steeringTorque_Nm = geom.numTires .* MzPerTire ./ ratioMag;
    handForce_N       = steeringTorque_Nm ./ geom.wheelRadius;

    out = struct('tireAngleDeg', tireAngleDeg, ...
                  'steeringWheelAngleDeg', steeringWheelAngleDeg, ...
                  'steeringRatio', ratio, ...
                  'MzPerTire_Nm', MzPerTire, ...
                  'steeringTorque_Nm', steeringTorque_Nm, ...
                  'handForce_N', handForce_N);

    thetaSW_atMax = steeringGeometryAngle(deg2rad(geom.maxTireRotationDeg), geom);
    predictedMaxWheelDeg = rad2deg(thetaSW_atMax);
    pctDiff = 100*(predictedMaxWheelDeg - geom.maxWheelRotationDeg)/geom.maxWheelRotationDeg;
    out.check = struct('predictedWheelAngleAtMaxTire_deg', predictedMaxWheelDeg, ...
                        'specifiedMaxWheelRotationDeg', geom.maxWheelRotationDeg, ...
                        'percentDifference', pctDiff);
end


function [thetaSW_rad, ratio] = steeringGeometryAngle(deltaTire_rad, geom)
%STEERINGGEOMETRYANGLE Tire steer angle -> steering-wheel angle + local ratio.
%   Idealized steering-arm crank: the arm (length armLength) is assumed
%   perpendicular to the tie rod at center, so rack travel
%   x(delta) = armLength*sin(delta). Pinion rotation = x*(2*pi/rackRatio)
%   [rackRatio in m/rev]; wheel/column angle = pinion angle*bevelGearRatio.
%   RATIO = d(thetaSW)/d(deltaTire) [rad/rad] via central finite difference
%   -- not constant, since this crank relationship is nonlinear. Replace
%   the one x_rack line below with your real CAD/measured toe curve if
%   you have one; nothing else needs to change.

    x_rack       = geom.armLength .* sin(deltaTire_rad);
    theta_pinion = x_rack .* (2*pi ./ geom.rackRatio);
    thetaSW_rad  = theta_pinion .* geom.bevelGearRatio;

    if nargout > 1
        h  = 1e-5; % rad
        f  = @(d) (geom.armLength.*sin(d)) .* (2*pi./geom.rackRatio) .* geom.bevelGearRatio;
        ratio = (f(deltaTire_rad+h) - f(deltaTire_rad-h)) / (2*h);
    end
end


function [Fy0, Mz0, info] = pacejkaPureLateral(alpha, Fz, gamma, P)
%PACEJKAPURELATERAL Pure-slip PAC2002/MF5.2 lateral force & aligning moment.
%   [FY0,MZ0] = PACEJKAPURELATERAL(ALPHA_rad, FZ_N, GAMMA_rad, P) evaluates
%   the standard published PAC2002 pure-slip equations (Pacejka, "Tyre and
%   Vehicle Dynamics"; MSC Adams/Tire "Using the PAC2002 Tire Model"). No
%   turn-slip / pressure terms (not present in a base PAC2002 file). Sign
%   convention follows whatever the .tir file was fitted in; for steering
%   effort you generally want abs(Mz0), which the caller applies.

    alpha = alpha(:);
    n = numel(alpha);
    if isscalar(Fz),    Fz    = repmat(Fz, n, 1);    else, Fz = Fz(:);       end
    if isscalar(gamma), gamma = repmat(gamma, n, 1); else, gamma = gamma(:); end

    R0   = P.UNLOADED_RADIUS;
    Fz0p = P.LFZ0 * P.FNOMIN;
    dfz  = (Fz - Fz0p) ./ Fz0p;

    % ---- Fy0 ----
    gamma_y = gamma .* P.LGAY;
    Cy  = P.PCY1 * P.LCY;
    muy = (P.PDY1 + P.PDY2.*dfz) .* (1 - P.PDY3.*gamma_y.^2) .* P.LMUY;
    Dy  = muy .* Fz;
    Ky  = P.PKY1.*Fz0p.*sin(2*atan(Fz./(P.PKY2*Fz0p))) .* (1-P.PKY3.*abs(gamma_y)) .* P.LFZ0 .* P.LKY;
    By  = Ky ./ safeDenom(Cy.*Dy);
    Shy = (P.PHY1+P.PHY2.*dfz).*P.LHY + P.PHY3.*gamma_y;
    Svy = Fz .* ((P.PVY1+P.PVY2.*dfz).*P.LVY + (P.PVY3+P.PVY4.*dfz).*gamma_y) .* P.LMUY;
    alpha_y = alpha + Shy;
    Ey = min((P.PEY1+P.PEY2.*dfz).*(1-(P.PEY3+P.PEY4.*gamma_y).*sign(alpha_y)).*P.LEY, 1);
    Bxarg = By.*alpha_y;
    Fy0 = Dy.*sin(Cy.*atan(Bxarg - Ey.*(Bxarg-atan(Bxarg)))) + Svy;

    % ---- Mz0 ----
    gamma_z = gamma .* P.LGAZ;
    Sht = P.QHZ1 + P.QHZ2.*dfz + (P.QHZ3+P.QHZ4.*dfz).*gamma_z;
    alpha_t = alpha + Sht;
    Shf = Shy + Svy ./ safeDenom(Ky);
    alpha_r = alpha + Shf;

    Bt = (P.QBZ1+P.QBZ2.*dfz+P.QBZ3.*dfz.^2) .* (1+P.QBZ4.*gamma_z+P.QBZ5.*abs(gamma_y)) .* (P.LKY/P.LMUY);
    Ct = P.QCZ1;
    Dt = Fz.*(P.QDZ1+P.QDZ2.*dfz).*(1+P.QDZ3.*gamma_z+P.QDZ4.*gamma_z.^2).*(R0/Fz0p).*P.LTR;
    Et = min((P.QEZ1+P.QEZ2.*dfz+P.QEZ3.*dfz.^2).*(1+(P.QEZ4+P.QEZ5.*gamma_z).*(2/pi).*atan(Bt.*Ct.*alpha_t)), 1);
    tArg = Bt.*alpha_t;
    t0 = Dt.*cos(Ct.*atan(tArg-Et.*(tArg-atan(tArg)))).*cos(alpha);

    Br = P.QBZ9*(P.LKY/P.LMUY) + P.QBZ10.*By.*Cy;
    Cr = 1; % fixed constant in the published formula (not QCZ1)
    Dr = Fz.*((P.QDZ6+P.QDZ7.*dfz).*P.LRES + (P.QDZ8+P.QDZ9.*dfz).*gamma_z).*R0.*P.LMUY;
    Mzr = Dr.*cos(Cr.*atan(Br.*alpha_r)).*cos(alpha);

    Mz0 = -t0.*Fy0 + Mzr;
    info = struct('By',By,'Cy',Cy,'Dy',Dy,'Ky',Ky,'Bt',Bt,'Dt',Dt,'Et',Et,'Br',Br,'Dr',Dr,'t0',t0,'Mzr',Mzr);
end

function d = safeDenom(d)
    % Avoid divide-by-zero without materially changing well-posed results.
    tol = 1e-6;
    idx = abs(d) < tol;
    s = sign(d(idx));
    s(s==0) = 1;
    d(idx) = tol .* s;
end


function P = readTIR(tirFile)
%READTIR Parse a PAC2002/MF5.2 .tir file into a struct (P.PCY1, P.FNOMIN, ...).
%   Ignores '!' and inline '$' comments and [SECTION] headers; every plain
%   "NAME = VALUE" line becomes a field (numeric, or string if quoted).

    if ~isfile(tirFile)
        error('readTIR:fileNotFound', 'Tire file not found: %s', tirFile);
    end
    fid = fopen(tirFile, 'r');
    if fid < 0
        error('readTIR:openFailed', 'Could not open tire file: %s', tirFile);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    P = struct();
    while true
        rawLine = fgetl(fid);
        if ~ischar(rawLine), break; end

        line = rawLine;
        bangIdx = strfind(line, '!');
        if ~isempty(bangIdx), line = line(1:bangIdx(1)-1); end
        dollarIdx = strfind(line, '$');
        if ~isempty(dollarIdx), line = line(1:dollarIdx(1)-1); end
        line = strtrim(line);

        if isempty(line) || startsWith(line, '[')
            continue
        end
        eqIdx = strfind(line, '=');
        if isempty(eqIdx)
            continue
        end
        eqIdx = eqIdx(1);

        name  = strtrim(line(1:eqIdx-1));
        value = strtrim(line(eqIdx+1:end));
        if isempty(name)
            continue
        end
        name = upper(name);
        if ~isvarname(name)
            continue
        end

        if startsWith(value, "'")
            P.(name) = strtrim(erase(value, "'"));
            continue
        end
        numVal = str2double(value);
        if isnan(numVal)
            continue
        end
        P.(name) = numVal;
    end
end
