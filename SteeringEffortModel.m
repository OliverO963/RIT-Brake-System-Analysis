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
%   4. Tire normal load (Fz) is driven by the SAME static + longitudinal +
%      lateral load-transfer formulas as the VD team's
%      Steering_Torque_Full_Calc.m (front axle: Fz_Fs = m*g*w_F/2 -
%      dFz_x/2, dFl = m*lat_g*g*h_cg*rd_F/track_F, Fz_OF/IF = Fz_Fs +/-
%      dFl), evaluated at that script's same 5 load cases. We do NOT
%      reproduce that script's combined-slip / torque-vectoring /
%      iterative steer-angle solve -- this model has no vehicle-speed or
%      cornering-equilibrium state (see point 1), so Fz is computed
%      per-case and then held fixed across the tire-angle sweep, same as
%      how a flat-trac Mz curve is measured at a fixed, specified load.
%   Not included: Ackermann angle difference, overturning moment (Mx),
%   scrub radius / kingpin offset, combined slip, torque vectoring.
%
% NEW IN THIS VERSION
%   - Fz now comes from 5 load cases (matching Steering_Torque_Full_Calc.m)
%     instead of a single nominal tire load.
%   - Sensitivity sweeps for 3 Brakes & Driver Controls parameters: total
%     wheel angular travel, intermediate column bevel-gear ratio, and
%     steering-wheel rim radius. Each gets its own figure (family-of-curves
%     vs wheel angle, plus peak torque/hand-force vs the swept parameter),
%     evaluated at the single load case that produces the highest baseline
%     peak torque (i.e. the worst case for driver effort).
%
% >>> EDIT THE "USER INPUT" SECTION BELOW <<<

clear; clc; close all;
set(0, 'DefaultFigureWindowStyle', 'docked');

%% ------------------------- USER INPUT -----------------------------
tirFile = "G:\Shared drives\RIT Formula SAE\Knowledge Center\Vehicle Dynamics\In-House VD Tools\Car Goals Models\Tire Models\AgileTireR20_AdjustedFX.tir";  % path to your PAC2002/MF5.2 .tir file

% VD-driven parameters
geom.armLength           = 0.0853;   % [m]     steering arm length (kingpin to tie-rod ball joint)
geom.rackTravel          = 0.0407; % [m] steering rack travel (0 degrees to full lock in one direction)
geom.maxTireRotationDeg  = (29.963 + 29.484)/2;      % [deg]   max tire steer angle from center, one direction (average of inner and outer to account for anti-Ackermann

% BDC-driven parameters (BASELINE values -- these are also the sweep
% baselines below: each sensitivity sweep varies ONE of these three while
% holding the other two at these baseline values)
geom.maxWheelRotationDeg = 100;     % [deg]   max steering-wheel rotation from center, one direction
geom.bevelGearRatio      = 1.0;    % [-]     column rotation / pinion rotation (1 = no bevel box)
geom.wheelRadius         = 0.115;    % [m]     steering wheel rim radius
geom.rackRatio           = geom.rackTravel/(geom.maxWheelRotationDeg/360);    % [m/rev] rack travel per ONE FULL revolution of the pinion

geom.Fz        = [];  % [N]   fallback only (LFZ0*FNOMIN); overwritten per load case below
geom.camberDeg = -1.5;   % [deg] static camber angle
geom.numTires  = 2;   % number of steered tires reacting torque through this shaft

nPoints = 200;        % resolution of the sweep from 0 to maxTireRotationDeg

% ---- Vehicle / normal-load (Fz) parameters -----------------------------
% Same values & formulas as Steering_Torque_Full_Calc.m's front-axle
% load-transfer block (m_total, w_front, h_cg_mm, wheelbase_mm, track_F_mm,
% K_roll_F/R_NmDeg). Only the load-transfer piece is reused -- not the
% combined-slip tire solve, TV algorithm, or iterative steer-angle solve.
vehicle.m_total     = 186 + 68;      % [kg]
vehicle.w_front     = 0.505;         % [-]  front weight distribution
vehicle.h_cg        = 271.5/1000;    % [m]
vehicle.wheelbase   = 1574.8/1000;   % [m]
vehicle.track_F     = 1219.2/1000;   % [m]
K_roll_F_NmDeg      = 971.1;         % [Nm/deg]
K_roll_R_NmDeg      = 971.1;         % [Nm/deg]
vehicle.roll_dist_F = K_roll_F_NmDeg/(K_roll_F_NmDeg + K_roll_R_NmDeg);

% How to collapse the outer/inner front tire loads into the single Fz this
% model applies to BOTH steered tires (this model, unlike the VD team's,
% does not solve outer/inner tires separately):
%   'average' -> Fz = mean(Fz_OF,Fz_IF). NOTE: lateral transfer cancels
%                EXACTLY in this mean (Fz_avg = Fz_Fs regardless of
%                lat_g), so cases sharing the same long_g collapse to the
%                same Fz. Longitudinal transfer still comes through.
%   'outer'   -> Fz = Fz_OF (the more heavily loaded, torque-dominant
%                tire). Captures the lateral component too, at the cost
%                of the "single averaged Fz" simplicity.
fzRepresentativeMode = 'average';

% Same 5 load cases as Steering_Torque_Full_Calc.m (name, lat_g, long_g)
loadCases = { ...
  'Skidpad steady',      1.00,  0.00; ...
  'Hairpin exit',        0.80,  0.80; ...
  'Hairpin exit hard',   0.60,  1.10; ...
  'Entry on regen',      0.80, -1.00; ...
  'Sweeper + regen',     1.00, -0.60};

% ---- Brakes & Driver Controls sensitivity sweeps -----------------------
wheelRotationSweepDeg  = linspace(90, 110, 5);   % [deg] total wheel travel from center, one direction
bevelGearRatioSweep    = linspace(0.8, 1.2, 5);  % [-]   +/-20% around baseline column/pinion ratio
wheelRadiusSweep_mm    = linspace(100, 130, 5);  % [mm]  steering wheel rim radius
%% ---------------------------------------------------------------------

P = readTIR(tirFile);
if isempty(geom.Fz)
    geom.Fz = P.LFZ0 * P.FNOMIN;
end

tireAngleSweepDeg = linspace(0.05, geom.maxTireRotationDeg, nPoints);

%% ----------------- Geometry consistency check (Fz-independent) --------
checkResult = steeringShaftTorque(tireAngleSweepDeg, P, geom);
fprintf('--- Geometry consistency check ---\n');
fprintf('At max tire angle (%.1f deg), this geometry gives a wheel angle of %.1f deg ', ...
        geom.maxTireRotationDeg, checkResult.check.predictedWheelAngleAtMaxTire_deg);
fprintf('(you specified %.1f deg, %.1f%% difference).\n', ...
        checkResult.check.specifiedMaxWheelRotationDeg, checkResult.check.percentDifference);
fprintf('If this is off, adjust armLength / rackRatio / bevelGearRatio.\n\n');

%% ----------------- Front-axle Fz per load case -------------------------
nCases = size(loadCases, 1);
caseFz = zeros(nCases, 3); % columns: Fz_OF, Fz_IF, Fz_avg
fprintf('--- Front-axle Fz by load case (load-transfer only, no TV/combined-slip) ---\n');
fprintf('%-19s %6s %6s %8s %8s %8s\n', 'case', 'lat_g', 'long_g', 'Fz_OF', 'Fz_IF', 'Fz_used');
for ii = 1:nCases
    [Fz_OF, Fz_IF, Fz_avg] = frontAxleFz(loadCases{ii,2}, loadCases{ii,3}, vehicle);
    caseFz(ii,:) = [Fz_OF, Fz_IF, Fz_avg];
    Fz_used = pickFz(caseFz(ii,:), fzRepresentativeMode);
    fprintf('%-19s %6.2f %6.2f %8.0f %8.0f %8.0f\n', ...
        loadCases{ii,1}, loadCases{ii,2}, loadCases{ii,3}, Fz_OF, Fz_IF, Fz_used);
end
if strcmp(fzRepresentativeMode, 'average')
    fprintf(['NOTE: fzRepresentativeMode = ''average'' -- lateral load transfer cancels\n' ...
             'in this mean, so cases with the same long_g show the same Fz_used above.\n' ...
             'Set fzRepresentativeMode = ''outer'' to keep the lateral component.\n']);
end
fprintf('\n');

%% ----------------- Baseline sweep, one curve per load case -------------
resultsByCase = cell(nCases, 1);
fprintf('--- Peak steering-shaft torque by load case (baseline BDC geometry) ---\n');
for ii = 1:nCases
    geomCase = geom;
    geomCase.Fz = pickFz(caseFz(ii,:), fzRepresentativeMode);
    resultsByCase{ii} = steeringShaftTorque(tireAngleSweepDeg, P, geomCase);
    [pk, idxPk] = max(resultsByCase{ii}.steeringTorque_Nm);
    fprintf('%-19s peak %6.2f N*m at tire angle %5.1f deg (wheel %5.1f deg), hand force %6.1f N\n', ...
        loadCases{ii,1}, pk, resultsByCase{ii}.tireAngleDeg(idxPk), ...
        resultsByCase{ii}.steeringWheelAngleDeg(idxPk), resultsByCase{ii}.handForce_N(idxPk));
end
fprintf('\n');

peakTorqueByCase = cellfun(@(r) max(r.steeringTorque_Nm), resultsByCase);
[worstPeakTorque, worstCaseIdx] = max(peakTorqueByCase);
worstCaseFz = pickFz(caseFz(worstCaseIdx,:), fzRepresentativeMode);
fprintf('Using "%s" (Fz = %.0f N) for the BDC sensitivity sweeps below -- it produces\n', ...
    loadCases{worstCaseIdx,1}, worstCaseFz);
fprintf('the highest baseline peak torque (%.2f N*m) of the 5 load cases.\n\n', worstPeakTorque);

%% ------------------------------- Plots (baseline, by load case) -------
caseLegend = loadCases(:,1);
co = lines(nCases);

figure('Name', 'Steering Effort by Load Case', 'Color', 'w', 'Position', [50 50 1000 700]);

subplot(2,2,1); hold on; grid on;
for ii = 1:nCases
    plot(resultsByCase{ii}.tireAngleDeg, resultsByCase{ii}.MzPerTire_Nm, 'Color', co(ii,:), 'LineWidth', 1.5);
end
xlabel('Tire steer angle [deg]'); ylabel('|M_z| per tire [N\cdotm]');
title('Tire self-aligning moment (Pacejka Mz_0)'); legend(caseLegend, 'Location', 'best');

subplot(2,2,2); hold on; grid on;
for ii = 1:nCases
    plot(resultsByCase{ii}.tireAngleDeg, resultsByCase{ii}.steeringWheelAngleDeg, 'Color', co(ii,:), 'LineWidth', 1.5);
end
xlabel('Tire steer angle [deg]'); ylabel('Steering wheel angle [deg]');
title('Steering system kinematics (Fz-independent -- lines overlap)'); legend(caseLegend, 'Location', 'best');

subplot(2,2,3); hold on; grid on;
for ii = 1:nCases
    plot(resultsByCase{ii}.steeringWheelAngleDeg, resultsByCase{ii}.steeringTorque_Nm, 'Color', co(ii,:), 'LineWidth', 1.5);
end
xlabel('Steering wheel angle [deg]'); ylabel('Steering shaft torque [N\cdotm]');
title('Torque driver must react, vs. wheel angle'); legend(caseLegend, 'Location', 'best');

subplot(2,2,4); hold on; grid on;
for ii = 1:nCases
    plot(resultsByCase{ii}.steeringWheelAngleDeg, resultsByCase{ii}.handForce_N, 'Color', co(ii,:), 'LineWidth', 1.5);
end
xlabel('Steering wheel angle [deg]'); ylabel('Equivalent rim hand force [N]');
title(sprintf('Hand force (rim radius = %.0f mm)', geom.wheelRadius*1000)); legend(caseLegend, 'Location', 'best');

sgtitle('Steering Effort by Load Case (baseline BDC geometry)');

%% ------------------- BDC sensitivity sweeps ----------------------------
sweepGeomBase = geom;
sweepGeomBase.Fz = worstCaseFz;

applyWheelRotation = @(g, v) setWheelRotation(g, v);
applyBevelRatio     = @(g, v) setBevelRatio(g, v);
applyWheelRadius    = @(g, v) setWheelRadius(g, v);

[wrResults, wrPeakT, wrPeakF] = runParamSweep(wheelRotationSweepDeg, applyWheelRotation, sweepGeomBase, P, tireAngleSweepDeg);
[bgResults, bgPeakT, bgPeakF] = runParamSweep(bevelGearRatioSweep, applyBevelRatio, sweepGeomBase, P, tireAngleSweepDeg);
[wrdResults, wrdPeakT, wrdPeakF] = runParamSweep(wheelRadiusSweep_mm, applyWheelRadius, sweepGeomBase, P, tireAngleSweepDeg);

plotSensitivity('Sensitivity to Total Wheel Angular Travel', wheelRotationSweepDeg, ...
    'Max wheel rotation, one direction [deg]', wrResults, wrPeakT, wrPeakF);

plotSensitivity('Sensitivity to Intermediate Column Bevel Gear Ratio', bevelGearRatioSweep, ...
    'Bevel gear ratio [-]', bgResults, bgPeakT, bgPeakF);

plotSensitivity('Sensitivity to Steering Wheel Radius', wheelRadiusSweep_mm, ...
    'Wheel rim radius [mm]', wrdResults, wrdPeakT, wrdPeakF);


%% ======================================================================
%  LOCAL FUNCTIONS (script-local functions, MATLAB R2016b+)
%  ======================================================================

function [Fz_OF, Fz_IF, Fz_avg] = frontAxleFz(lat_g, long_g, vehicle)
%FRONTAXLEFZ Front-axle outer/inner tire normal load from static +
%   longitudinal + lateral load transfer. Same formulas as the front-axle
%   block of Steering_Torque_Full_Calc.m's solve_corner (load transfer
%   only -- no combined-slip / TV / steer-angle iteration).
    g = 9.81;
    dFz_x = vehicle.m_total * long_g * g * vehicle.h_cg / vehicle.wheelbase;
    Fz_Fs = vehicle.m_total * g * vehicle.w_front/2 - dFz_x/2;
    dFl   = vehicle.m_total * lat_g * g * vehicle.h_cg * vehicle.roll_dist_F / vehicle.track_F;
    Fz_OF = max(Fz_Fs + dFl, 0);
    Fz_IF = max(Fz_Fs - dFl, 0);
    Fz_avg = (Fz_OF + Fz_IF) / 2;
end

function Fz = pickFz(caseFzRow, mode)
%PICKFZ Collapse [Fz_OF Fz_IF Fz_avg] to a single Fz per fzRepresentativeMode.
    switch mode
        case 'average', Fz = caseFzRow(3);
        case 'outer',   Fz = caseFzRow(1);
        otherwise, error('pickFz:badMode', 'fzRepresentativeMode must be ''average'' or ''outer''.');
    end
end

function g = setWheelRotation(g, maxWheelRotationDeg)
%SETWHEELROTATION Update total wheel travel and its dependent rack ratio.
    g.maxWheelRotationDeg = maxWheelRotationDeg;
    g.rackRatio = g.rackTravel / (maxWheelRotationDeg/360);
end

function g = setBevelRatio(g, bevelGearRatio)
    g.bevelGearRatio = bevelGearRatio;
end

function g = setWheelRadius(g, wheelRadius_mm)
    g.wheelRadius = wheelRadius_mm / 1000;
end

function [sweepResults, peakTorque_Nm, peakHandForce_N] = runParamSweep(paramValues, applyFcn, baseGeom, P, tireAngleSweepDeg)
%RUNPARAMSWEEP Run steeringShaftTorque once per value in PARAMVALUES,
%   modifying baseGeom via APPLYFCN(geom, value) each time.
    n = numel(paramValues);
    sweepResults = cell(n, 1);
    peakTorque_Nm = zeros(n, 1);
    peakHandForce_N = zeros(n, 1);
    for k = 1:n
        g = applyFcn(baseGeom, paramValues(k));
        r = steeringShaftTorque(tireAngleSweepDeg, P, g);
        sweepResults{k} = r;
        peakTorque_Nm(k) = max(r.steeringTorque_Nm);
        peakHandForce_N(k) = max(r.handForce_N);
    end
end

function plotSensitivity(figTitle, paramValues, paramLabel, sweepResults, peakTorque_Nm, peakHandForce_N)
%PLOTSENSITIVITY One figure per BDC parameter: family-of-curves (torque
%   and hand force vs wheel angle, one line per swept value) plus
%   peak-torque / peak-hand-force vs the swept parameter.
    n = numel(paramValues);
    co = lines(n);
    legendLabels = arrayfun(@(v) sprintf('%.3g', v), paramValues, 'UniformOutput', false);

    figure('Name', figTitle, 'Color', 'w', 'Position', [100 100 1000 700]);

    subplot(2,2,1); hold on; grid on;
    for k = 1:n
        plot(sweepResults{k}.steeringWheelAngleDeg, sweepResults{k}.steeringTorque_Nm, ...
            'Color', co(k,:), 'LineWidth', 1.5);
    end
    xlabel('Steering wheel angle [deg]'); ylabel('Steering shaft torque [N\cdotm]');
    title('Torque vs. wheel angle'); legend(legendLabels, 'Location', 'best');

    subplot(2,2,2); hold on; grid on;
    for k = 1:n
        plot(sweepResults{k}.steeringWheelAngleDeg, sweepResults{k}.handForce_N, ...
            'Color', co(k,:), 'LineWidth', 1.5);
    end
    xlabel('Steering wheel angle [deg]'); ylabel('Hand force [N]');
    title('Hand force vs. wheel angle'); legend(legendLabels, 'Location', 'best');

    subplot(2,2,3);
    plot(paramValues, peakTorque_Nm, '-o', 'LineWidth', 1.5); grid on;
    xlabel(paramLabel); ylabel('Peak steering torque [N\cdotm]');
    title('Peak torque sensitivity');

    subplot(2,2,4);
    plot(paramValues, peakHandForce_N, '-o', 'LineWidth', 1.5); grid on;
    xlabel(paramLabel); ylabel('Peak hand force [N]');
    title('Peak hand force sensitivity');

    sgtitle(figTitle);
end

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