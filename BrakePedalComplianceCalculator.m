%% RIT Racing Brake Pedal Compliance Calculator
% The objective of this script is to quantify the effects of brake pedal
% compliance on brake system performance.
% Created by Oliver Owen 22 April 2026
% Last updated by Oliver Owen on 02 June 2026

clear; clc; close all;
colors = ["#c44536", "#197278", "#772e25", "#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#a0c4ff", "#bdb2ff", "#ffc6ff"];

%% Global System Parameters
global frontMC_A rearMC_A frontPiston_A rearPiston_A frontPiston_num rearPiston_num bias pedalRatio lockPres
frontMC_BD = 0.7; % front master cylinder bore diameter [in]
rearMC_BD = 0.7; % rear master cylinder bore diameter [in]
frontMC_A = (pi/4)*(frontMC_BD^2); % front master cylinder bore area [in2]
rearMC_A = (pi/4)*(rearMC_BD^2); % rear master cylinder bore area [in2]

frontPiston_D = 0.619; % front piston diameter [in]
rearPiston_D = 0.558; % rear piston diameter [in]
frontPiston_A = (pi/4)*(frontPiston_D^2); % front piston area [in2]
rearPiston_A = (pi/4)*(rearPiston_D^2); % rear piston area [in2]

frontPiston_num = 6; % front piston count [-]
rearPiston_num = 4; % rear piston count [-]

bias = 0.6; % front brake bias at pedal [-]

lockForce = 200; % pedal input force to lock brakes [lbf]
pedalRatio = 6.001; % pedal ratio [-]
lockPres = lockForce*pedalRatio; % line pressure to lock

mcLength = 6.1528; % resting master cylinder length [in]
pedalLength = 8; % pedal length [in]

padMu = 0.5; % brake pad coefficient of friction [-]
padApplyD = 7.3; % diameter of circle through pad centerlne [in]

gasSpringForce = 6.18; % fixed gas spring reaction force [lbf], update to make variable with pedal angle

%% Error Bound Parameters (+/- 10% on uncertain deformation estimates)
errorFactor = 0.10; % fractional uncertainty on pedal, caliper, and bias bar max deformations [-]

%% Main Loop — Nominal, Upper Bound (+10%), Lower Bound (-10%)
forceVec = linspace(0, lockForce, 1000);
% forceVec = 200;

for n = 1:1:length(forceVec)
    if (forceVec(n) - gasSpringForce > 0)
        [frontHoopComp(n), rearHoopComp(n)] = hoopDef(forceVec(n) - gasSpringForce); % hard line hoop deformation contribution to master cylinder length change [in]
        [frontSoftHoopComp(n), rearSoftHoopComp(n)] = softHoopDef(forceVec(n) - gasSpringForce); % hoop deformation contribution to master cylinder length change [in]
        [frontPadTakeupComp(n), rearPadTakeupComp(n)] = padTakeup(); % pad takeup contribution to master cylinder length change [in]
        [frontPadCompComp(n), rearPadCompComp(n)] = padComp(forceVec(n)); % pad compression
        [frontFluidComp(n), rearFluidComp(n)] = fluidComp(forceVec(n)); % fluid compression
    
        % --- Nominal caliper/bias bar deformation ---
        [frontPartsComp(n), rearPartsComp(n)] = parts(forceVec(n) - gasSpringForce, 1.0); % caliper and bias bar deformation contribution to master cylinder length change [in]
        % --- Upper bound (+10%) caliper/bias bar deformation ---
        [frontPartsCompHigh(n), rearPartsCompHigh(n)] = parts(forceVec(n) - gasSpringForce, 1.0 + errorFactor);
        % --- Lower bound (-10%) caliper/bias bar deformation ---
        [frontPartsCompLow(n), rearPartsCompLow(n)] = parts(forceVec(n) - gasSpringForce, 1.0 - errorFactor);        
    
        frontMCFreePlayComp(n) = 0.04; % front master cylinder free play [in]
        rearMCFreePlayComp(n) = 0.047; % rear master cylinder free play [in]
        biasBarPlayComp(n) = 0.01; % bias bar free play [in]
    else
        % nothing if gas spring reacts input force
        frontHoopComp(n) = 0;
        rearHoopComp(n) = 0;

        frontSoftHoopComp(n) = 0;
        rearSoftHoopComp(n) = 0;

        frontPadTakeupComp(n) = 0;
        rearPadTakeupComp(n) = 0;

        frontPartsComp(n) = 0;
        rearPartsComp(n) = 0;
        frontPartsCompHigh(n) = 0;
        rearPartsCompHigh(n) = 0;
        frontPartsCompLow(n) = 0;
        rearPartsCompLow(n) = 0;

        frontPadCompComp(n) = 0; 
        rearPadCompComp(n) = 0; 

        frontFluidComp(n) = 0;
        rearFluidComp(n) = 0;

        frontMCFreePlayComp(n) = 0; 
        rearMCFreePlayComp(n) = 0; 

        biasBarPlayComp(n) = 0; 
    end

    % --- Nominal total compliance ---
    totalFrontComp(n) = frontHoopComp(n) + frontSoftHoopComp(n) + frontPadTakeupComp(n) + frontPadCompComp(n) + frontFluidComp(n) + frontPartsComp(n) + frontMCFreePlayComp(n) + biasBarPlayComp(n); % total front master cylinder length change [in]
    totalRearComp(n)  = rearHoopComp(n) + rearSoftHoopComp(n) + rearPadTakeupComp(n)  + rearPadCompComp(n) + rearFluidComp(n) + rearPartsComp(n)  + rearMCFreePlayComp(n)  + biasBarPlayComp(n);  % total rear master cylinder length change [in]

    % --- Upper bound total compliance (higher deformations → more MC travel) ---
    totalFrontCompHigh(n) = frontHoopComp(n) + frontSoftHoopComp(n) + frontPadTakeupComp(n) + frontPadCompComp(n) + frontFluidComp(n) + frontPartsCompHigh(n) + frontMCFreePlayComp(n) + biasBarPlayComp(n);
    totalRearCompHigh(n)  = rearHoopComp(n) + rearSoftHoopComp(n)  + rearPadTakeupComp(n)  + rearPadCompComp(n) + rearFluidComp(n)  + rearPartsCompHigh(n)  + rearMCFreePlayComp(n)  + biasBarPlayComp(n);

    % --- Lower bound total compliance ---
    totalFrontCompLow(n) = frontHoopComp(n) + frontSoftHoopComp(n) + frontPadTakeupComp(n) + frontPadCompComp(n) + frontFluidComp(n) + frontPartsCompLow(n) + frontMCFreePlayComp(n) + biasBarPlayComp(n);
    totalRearCompLow(n)  = rearHoopComp(n) + rearSoftHoopComp(n)  + rearPadTakeupComp(n)  + rearPadCompComp(n) + rearFluidComp(n)  + rearPartsCompLow(n)  + rearMCFreePlayComp(n)  + biasBarPlayComp(n);

    % --- Nominal pedal travel ---
    totalAverageComp(n) = (totalFrontComp(n) + totalRearComp(n)) / 2;
    newMClength(n) = mcLength - totalAverageComp(n);
    angPedalTravel(n) = 89.959 - rad2deg(acos((1.367^2+6^2-newMClength(n)^2)/(2*1.367*6)));

    pedalMaxDef = 0.03; % maximum pedal deformation at end [in] — nominal
    pedalDef(n) = (pedalMaxDef / lockForce) * (forceVec(n) - gasSpringForce);
    totalPedalTravel(n) = (pedalLength * sin(deg2rad(angPedalTravel(n)))) + pedalDef(n); % nominal pedal travel [in]

    % --- Upper bound pedal travel (+10% on pedal max def and parts deformation) ---
    pedalMaxDefHigh = pedalMaxDef * (1.0 + errorFactor);
    pedalDefHigh(n) = (pedalMaxDefHigh / lockForce) * (forceVec(n) - gasSpringForce);
    totalAverageCompHigh(n) = (totalFrontCompHigh(n) + totalRearCompHigh(n)) / 2;
    newMClengthHigh(n) = mcLength - totalAverageCompHigh(n);
    angPedalTravelHigh(n) = 89.959 - rad2deg(acos((1.367^2+6^2-newMClengthHigh(n)^2)/(2*1.367*6)));
    totalPedalTravelHigh(n) = (pedalLength * sin(deg2rad(angPedalTravelHigh(n)))) + pedalDefHigh(n); % upper bound pedal travel [in]

    % --- Lower bound pedal travel (-10% on pedal max def and parts deformation) ---
    pedalMaxDefLow = pedalMaxDef * (1.0 - errorFactor);
    pedalDefLow(n) = (pedalMaxDefLow / lockForce) * (forceVec(n) - gasSpringForce);
    totalAverageCompLow(n) = (totalFrontCompLow(n) + totalRearCompLow(n)) / 2;
    newMClengthLow(n) = mcLength - totalAverageCompLow(n);
    angPedalTravelLow(n) = 89.959 - rad2deg(acos((1.367^2+6^2-newMClengthLow(n)^2)/(2*1.367*6)));
    totalPedalTravelLow(n) = (pedalLength * sin(deg2rad(angPedalTravelLow(n)))) + pedalDefLow(n); % lower bound pedal travel [in]

    frontBrakeTorque(n) = ((forceVec(n)*pedalRatio*bias)/frontMC_A)*frontPiston_A*frontPiston_num*padMu*(padApplyD/2);
    rearBrakeTorque(n)  = ((forceVec(n)*pedalRatio*(1-bias))/rearMC_A)*rearPiston_A*rearPiston_num*padMu*(padApplyD/2);
end

figure(1);
hold on;
title("Linear Pedal Travel At End and Brake Torque vs. Input Force");
xlabel("Force [lbf]");

yyaxis left;
ylabel("Pedal Travel at End [in]");

% Shaded error band for pedal travel (fill between high and low bounds)
fillX = [forceVec, fliplr(forceVec)];
fillY = [totalPedalTravelHigh, fliplr(totalPedalTravelLow)];
fill(fillX, fillY, sscanf(colors(1),'#%2x%2x%2x',[1 3])/255, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', "\pm10% Metallic Part Deformation Uncertainty");

% Nominal pedal travel line
plot(forceVec, totalPedalTravel, 'Color', colors(1), 'LineWidth', 2, 'DisplayName', "Pedal Travel (Nominal)");

% Upper and lower bound lines (dashed)
plot(forceVec, totalPedalTravelHigh, 'Color', colors(1), 'LineWidth', 1, 'LineStyle', '--', 'DisplayName', "Pedal Travel (+10%)");
plot(forceVec, totalPedalTravelLow,  'Color', colors(1), 'LineWidth', 1, 'LineStyle', ':', 'DisplayName', "Pedal Travel (-10%)");

axis([0, lockForce, 0, inf]);
text(forceVec(end) - 40, totalPedalTravel(end) - 0.3, "Max Pedal Travel = " + totalPedalTravel(end) + " in", 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors(1));

% yyaxis right;
% ylabel("Braking Torque [lbf*in]");
% plot(forceVec, frontBrakeTorque, 'Color', colors(2), 'LineWidth', 2, 'DisplayName', "Front Brake Torque");
% text(forceVec(end) - 60, frontBrakeTorque(end) + 100, "Max Front Brake Torque = " + frontBrakeTorque(end) + " lbf*in", 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors(2));
% plot(forceVec, rearBrakeTorque, 'Color', colors(3), 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', "Rear Brake Torque");
% text(forceVec(end) - 60, rearBrakeTorque(end) + 100, "Max Rear Brake Torque = " + rearBrakeTorque(end) + " lbf*in", 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', colors(3));

%% Experimental Test Data (corrected for sine error of dial indicator)
testForce_t1   = [0, 5.82, 25, 50, 75, 100]; % Trial 1 scale force [lbf]
testDisp_t1    = [0, 0.119, 0.617, 0.763, 0.864, 0.961]; % Trial 1 corrected displacement [in]

testForce_t2   = [0, 6.6,  25, 50, 75, 100]; % Trial 2 scale force [lbf]
testDisp_t2    = [0, 0.119, 0.616, 0.761, 0.864, 0.951]; % Trial 2 corrected displacement [in]

testForce_t3   = [0, 7.45, 25, 50, 75, 100]; % Trial 3 scale force [lbf]
testDisp_t3    = [0, 0.125, 0.600, 0.747, 0.872, 0.950]; % Trial 3 corrected displacement [in]

yyaxis left;
scatter(testForce_t1, testDisp_t1, 40, sscanf(colors(5),'#%2x%2x%2x',[1 3])/255, 'o', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', "Test Data — Trial 1");
scatter(testForce_t2, testDisp_t2, 40, sscanf(colors(6),'#%2x%2x%2x',[1 3])/255, 's', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', "Test Data — Trial 2");
scatter(testForce_t3, testDisp_t3, 40, sscanf(colors(7),'#%2x%2x%2x',[1 3])/255, '^', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', "Test Data — Trial 3");

axes = gca;
axes.YAxis(1).Color = "k";
axes.YAxis(2).Color = "k";
grid on;
legend('Location', 'northwest');
fontname("Montserrat");

%% Debug Output — Max Values at lockForce (200 lbf)
fprintf('\n========== DEBUG: Max Component Values at %.0f lbf ==========\n', lockForce);
fprintf('%-45s %10s %10s\n', 'Component', 'Front [in]', 'Rear [in]');
fprintf('%s\n', repmat('-', 1, 67));
fprintf('%-45s %10.6f %10.6f\n', 'Hoop Deformation',           frontHoopComp(end),        rearHoopComp(end));
fprintf('%-45s %10.6f %10.6f\n', 'Pad Takeup',                 frontPadTakeupComp(end),   rearPadTakeupComp(end));
fprintf('%-45s %10.6f %10.6f\n', 'Pad Compression',    frontPadCompComp(end),     rearPadCompComp(end));
fprintf('%-45s %10.6f %10.6f\n', 'Fluid Compression',    frontFluidComp(end),     rearFluidComp(end));
fprintf('%-45s %10.6f %10.6f\n', 'Caliper + Bias Bar Def',     frontPartsComp(end),       rearPartsComp(end));
fprintf('%-45s %10.6f %10.6f\n', 'MC Free Play',               frontMCFreePlayComp(end),  rearMCFreePlayComp(end));
fprintf('%-45s %10.6f %10.6f\n', 'Bias Bar Free Play',         biasBarPlayComp(end),      biasBarPlayComp(end));
fprintf('%s\n', repmat('-', 1, 67));
fprintf('%-45s %10.6f %10.6f\n', 'Total Compliance',           totalFrontComp(end),       totalRearComp(end));
fprintf('%-45s %10.6f\n',        'Average Total Compliance',   totalAverageComp(end));
fprintf('%s\n', repmat('-', 1, 67));
fprintf('%-45s %10.6f\n',        'New MC Length [in]',         newMClength(end));
fprintf('%-45s %10.6f\n',        'Angular Pedal Travel [deg]', angPedalTravel(end));
fprintf('%-45s %10.6f\n',        'Pedal Structural Deflection [in]', pedalDef(end));
fprintf('%-45s %10.6f\n',        'Total Pedal Travel [in]',    totalPedalTravel(end));
fprintf('%-45s %10.6f\n',        'Total Pedal Travel +10%% [in]', totalPedalTravelHigh(end));
fprintf('%-45s %10.6f\n',        'Total Pedal Travel -10%% [in]', totalPedalTravelLow(end));
fprintf('=============================================================\n\n');

%% Debug Output — Max Values at 0 lbf
fprintf('\n========== DEBUG: Max Component Values at %.0f lbf ==========\n', 0);
fprintf('%-45s %10s %10s\n', 'Component', 'Front [in]', 'Rear [in]');
fprintf('%s\n', repmat('-', 1, 67));
fprintf('%-45s %10.6f %10.6f\n', 'Hoop Deformation',           frontHoopComp(1),        rearHoopComp(1));
fprintf('%-45s %10.6f %10.6f\n', 'Pad Takeup',                 frontPadTakeupComp(1),   rearPadTakeupComp(1));
fprintf('%-45s %10.6f %10.6f\n', 'Pad Compression',    frontPadCompComp(1),     rearPadCompComp(1));
fprintf('%-45s %10.6f %10.6f\n', 'Fluid Compression',    frontFluidComp(1),     rearFluidComp(1));
fprintf('%-45s %10.6f %10.6f\n', 'Caliper + Bias Bar Def',     frontPartsComp(1),       rearPartsComp(1));
fprintf('%-45s %10.6f %10.6f\n', 'MC Free Play',               frontMCFreePlayComp(1),  rearMCFreePlayComp(1));
fprintf('%-45s %10.6f %10.6f\n', 'Bias Bar Free Play',         biasBarPlayComp(1),      biasBarPlayComp(1));
fprintf('%s\n', repmat('-', 1, 67));
fprintf('%-45s %10.6f %10.6f\n', 'Total Compliance',           totalFrontComp(1),       totalRearComp(1));
fprintf('%-45s %10.6f\n',        'Average Total Compliance',   totalAverageComp(1));
fprintf('%s\n', repmat('-', 1, 67));
fprintf('%-45s %10.6f\n',        'New MC Length [in]',         newMClength(1));
fprintf('%-45s %10.6f\n',        'Angular Pedal Travel [deg]', angPedalTravel(1));
fprintf('%-45s %10.6f\n',        'Pedal Structural Deflection [in]', pedalDef(1));
fprintf('%-45s %10.6f\n',        'Total Pedal Travel [in]',    totalPedalTravel(1));
fprintf('%-45s %10.6f\n',        'Total Pedal Travel +10%% [in]', totalPedalTravelHigh(1));
fprintf('%-45s %10.6f\n',        'Total Pedal Travel -10%% [in]', totalPedalTravelLow(1));
fprintf('=============================================================\n\n');

%% Hardline Hoop Deformation Function
function [frontMC_delta, rearMC_delta] = hoopDef(force)
    % Caculates change in master cylinder length due to hardline hoop deformation.
    % Input: pedal force [lbf]
    % Outputs: front master cylinder linear displacement [in], rear master cylinder linear displacement [in]

    % pull globals
    global frontMC_A rearMC_A pedalRatio bias

    % internal function inputs, hard line paramters
    line_OD = 0.1875; % outer diameter [in]
    line_t = 0.028; % wall thickness [in]
    line_fL = 47.8; % front length [in]
    line_rL = 112.5; % rear length [in]
    line_nu = 0.29; % Poisson's ratio [in]
    line_E = 28000000; % Young's modulus [psi]
    
    line_ID = line_OD - 2*line_t; % line inner diameter [in]
    line_r = line_ID/2; % line inner radius [in]

    frontPres = (force*pedalRatio*bias)/frontMC_A; % front system line pressure [psi]
    rearPres = (force*pedalRatio*(1-bias))/rearMC_A; % rear system line pressure [psi]

    frontHoopStress = (frontPres*line_ID)/(2*line_t); % front line hoop stress [psi]
    frontAxialStress = frontHoopStress/2; % front line axial stress [psi];
    rearHoopStress = (rearPres*line_ID)/(2*line_t); % rear line hoop stress [psi]
    rearAxialStress = rearHoopStress/2; % rear line axial stress [psi];

    frontStrain = (1/line_E)*(frontHoopStress - (line_nu*frontAxialStress)); % front line strain [in/in]
    rearStrain = (1/line_E)*(rearHoopStress - (line_nu*rearAxialStress)); % rear line strain [in/in]

    frontDefRad = line_r + line_r*frontStrain; % front deformed radius [in]
    rearDefRad = line_r + line_r*rearStrain; % rear deformed radius [in]

    frontVol = pi*line_fL*(frontDefRad^2 - line_r^2); % front volume demand due to hoop strain [in3]
    rearVol = pi*line_rL*(rearDefRad^2 - line_r^2); % rear volume demand due to hoop strain [in3]

    frontMC_delta = frontVol/frontMC_A; % front master cylinder change in length due to hoop strain [in]
    rearMC_delta = rearVol/rearMC_A; % rear master cylinder change in length due to hoop strain [in]
end

%% Softline Hoop Deformation Function
function [frontMC_delta, rearMC_delta] = softHoopDef(force)
    % Caculates change in master cylinder length due to hardline hoop deformation.
    % Input: pedal force [lbf]
    % Outputs: front master cylinder linear displacement [in], rear master cylinder linear displacement [in]

    % pull globals
    global frontMC_A rearMC_A pedalRatio bias

    % internal function inputs, hard line paramters
    line_OD = 0.3; % outer diameter [in]
    line_t = 0.0875; % wall thickness [in]
    line_fL = 64; % front length [in]
    line_rL = 42; % rear length [in]
    line_nu = 0.35; % Poisson's ratio [in]
    line_E = 3000000; % Young's modulus [psi]
    
    line_ID = line_OD - 2*line_t; % line inner diameter [in]
    line_r = line_ID/2; % line inner radius [in]

    frontPres = (force*pedalRatio*bias)/frontMC_A; % front system line pressure [psi]
    rearPres = (force*pedalRatio*(1-bias))/rearMC_A; % rear system line pressure [psi]

    frontHoopStress = (frontPres*line_ID)/(2*line_t); % front line hoop stress [psi]
    frontAxialStress = frontHoopStress/2; % front line axial stress [psi];
    rearHoopStress = (rearPres*line_ID)/(2*line_t); % rear line hoop stress [psi]
    rearAxialStress = rearHoopStress/2; % rear line axial stress [psi];

    frontStrain = (1/line_E)*(frontHoopStress - (line_nu*frontAxialStress)); % front line strain [in/in]
    rearStrain = (1/line_E)*(rearHoopStress - (line_nu*rearAxialStress)); % rear line strain [in/in]

    frontDefRad = line_r + line_r*frontStrain; % front deformed radius [in]
    rearDefRad = line_r + line_r*rearStrain; % rear deformed radius [in]

    frontVol = pi*line_fL*(frontDefRad^2 - line_r^2); % front volume demand due to hoop strain [in3]
    rearVol = pi*line_rL*(rearDefRad^2 - line_r^2); % rear volume demand due to hoop strain [in3]

    frontMC_delta = frontVol/frontMC_A; % front master cylinder change in length due to hoop strain [in]
    rearMC_delta = rearVol/rearMC_A; % rear master cylinder change in length due to hoop strain [in]
end


%% Pad Takeup Function
function [frontMC_delta, rearMC_delta] = padTakeup()
    % Calculates change in master cylinder length due to pad takeup.
    % Inputs: none (all globals)
    % Outputs: front master cylinder linear displacement [in], rear master cylinder linear displacement [in]

    global frontMC_A rearMC_A frontPiston_A rearPiston_A frontPiston_num rearPiston_num
    
    % internal function inputs
    frontPiston_travel = 0.015; % travel of a single front piston to meet the pad [in]
    rearPiston_travel = 0.002; % travel of a single rear piston to meet the pad [in]

    frontVol = frontPiston_travel*frontPiston_A*frontPiston_num; % front volume demand due to pad takeup [in3]
    rearVol = rearPiston_travel*rearPiston_A*rearPiston_num; % rear volume demand due to pad takeup [in3]

    frontMC_delta = frontVol/frontMC_A; % front master cylinder change in length due to pad takeup [in]
    rearMC_delta = rearVol/rearMC_A; % rear master cylinder change in length due to pad takeup [in]
end

%% Brake Fluid Compression 
function [frontMC_delta, rearMC_delta] = fluidComp(force)
    % Calculates change in master cylinder length due to brake fluid compressibility.
    % Inputs: force
    % Outputs: front master cylinder linear displacement [in], rear master cylinder linear displacement [in]

    global frontMC_A rearMC_A pedalRatio bias

    airFactor = 1;

    % Hard line geometry
    hLine_OD = 0.1875;  % outer diameter [in]
    hLine_t  = 0.028;   % wall thickness [in]
    hLine_fL = 47.8;    % front hard line length [in]
    hLine_rL = 112.5;   % rear hard line length [in]
    hLine_ID = hLine_OD - 2*hLine_t;

    % Soft line geometry
    sLine_OD = 0.3;        % soft line outer diameter [in]
    sLine_t  = 0.0875;        % soft line wall thickness [in]
    sLine_fL = 64;        % front soft line length [in]
    sLine_rL = 42;        % rear soft line length [in]
    sLine_ID = sLine_OD - 2*sLine_t;

    % Caliper internal volume - not included for now
    frontCaliperVol = 0;  % front caliper internal fluid volume [in3]
    rearCaliperVol  = 0;  % rear caliper internal fluid volume [in3]

    % Total fluid volume
    frontLine_V = (pi/4) * hLine_ID^2 * hLine_fL ...
                + (pi/4) * sLine_ID^2 * sLine_fL ...
                + frontCaliperVol;  % total front circuit fluid volume [in3]
    rearLine_V  = (pi/4) * hLine_ID^2 * hLine_rL ...
                + (pi/4) * sLine_ID^2 * sLine_rL ...
                + rearCaliperVol;   % total rear circuit fluid volume [in3]

    B = 19535 * 14.5038;  % bulk modulus [bar] -> [psi]
    B = B * airFactor;

    frontPres = (force * pedalRatio * bias)       / frontMC_A;  % front line pressure [psi]
    rearPres  = (force * pedalRatio * (1 - bias)) / rearMC_A;   % rear line pressure [psi]

    % Volumetric compression: dV = V * dP / B
    frontVol_delta = frontLine_V * (frontPres / B);  % front volume change [in3]
    rearVol_delta  = rearLine_V  * (rearPres  / B);  % rear volume change [in3]

    frontMC_delta = frontVol_delta / frontMC_A;  % front MC length change [in]
    rearMC_delta  = rearVol_delta  / rearMC_A;   % rear MC length change [in]
end
%% Pad Compression Function
function [frontMC_delta, rearMC_delta] = padComp(force)
    % Calculates change in master cylinder length due to pad compression.
    % Inputs: pedal force [lbf]
    % Outputs: front master cylinder linear displacement [in], rear master cylinder linear displacement [in]

    global frontMC_A rearMC_A frontPiston_A rearPiston_A frontPiston_num rearPiston_num pedalRatio bias

    % internal function inputs
    frontPad_A = 2.36; % front pad area [in2]
    rearPad_A = 1.015; % rear pad area [in2]
    padThick = 0.3; % pad thickness [in]
    pad_Ec = 50000; % pad material compressive modulus [psi]

    frontPres = (force*pedalRatio*bias)/frontMC_A; % front system line pressure [psi]
    rearPres = (force*pedalRatio*(1-bias))/rearMC_A; % rear system line pressure [psi]

    frontCompForce = frontPres*(frontPiston_A*(frontPiston_num/2)); % single front pad compressing force [lbf]
    rearCompForce = rearPres*(rearPiston_A*(rearPiston_num/2)); % single rear pad compressing force [lbf]

    frontPad_Stress = frontCompForce/frontPad_A; % single front pad stress [psi]
    rearPad_Stress = rearCompForce/rearPad_A; % single rear pad stress [psi]

    frontPad_Strain = frontPad_Stress/pad_Ec; % single front pad strain [in/in]
    rearPad_Strain = rearPad_Stress/pad_Ec; % single rear pad strain [in/in]

    frontPad_newThick = padThick*(1-frontPad_Strain); % new front pad thickness [in]
    rearPad_newThick = padThick*(1-rearPad_Strain); % new front pad thickness [in]

    frontPad_delta = padThick - frontPad_newThick; % change in front pad thickness [in]
    rearPad_delta = padThick - rearPad_newThick; % change in rear pad thickness [in]

    frontVol = frontPiston_A*frontPiston_num*frontPad_delta; % front volume demand [in]
    rearVol = rearPiston_A*rearPiston_num*rearPad_delta; % rear volume demand [in]

    frontMC_delta = frontVol/frontMC_A; % front master cylinder change in length due to pad compression [in]
    rearMC_delta = rearVol/rearMC_A; % rear master cylinder change in length due to pad compression [in]
end

%% Caliper and Bias Bar Deformation Function
function [frontMC_delta, rearMC_delta] = parts(force, defScale)
    % Calculates change in master cylinder length due to caliper and bias bar deformation
    % Inputs: pedal force [lbf], defScale [-] (1.0 = nominal, 1.1 = +10%, 0.9 = -10%)
    % Outputs: front master cylinder linear displacement [in], rear master cylinder linear displacement [in]

    global frontMC_A rearMC_A frontPiston_A rearPiston_A frontPiston_num rearPiston_num pedalRatio bias lockPres

    % internal function inputs (scaled by defScale for uncertainty bounding)
    frontCaliper_maxDef = 0.0015 * defScale; % front caliper maximum deformation [in]
    rearCaliper_maxDef  = 0.0015 * defScale; % rear caliper maximum deformation [in]
    biasBar_maxDef      = 0.008  * defScale; % bias bar max deformation [in]
    mcMount_maxDef      = 0.033 * defScale;   % master cylinder pivot max deformation [in]

    frontPres = (force*pedalRatio*bias)/frontMC_A; % front system line pressure [psi]
    rearPres = (force*pedalRatio*(1-bias))/rearMC_A; % rear system line pressure [psi]

    frontCaliper_def = (frontCaliper_maxDef/lockPres) * frontPres; % front caliper deformation due to current applied force [in]
    rearCaliper_def = (rearCaliper_maxDef/lockPres) * rearPres; % rear caliper deformation due to current applied force [in]
    biasBar_def = (biasBar_maxDef/lockPres) * (frontPres+rearPres); % bias bar deformation due to current applied force [in]
    mcMount_def = (mcMount_maxDef/lockPres) * (frontPres+rearPres); % master cylinder mount deformation due to current applied force [in]

    frontCaliperVol = frontCaliper_def*frontPiston_A*frontPiston_num; % front caliper volume demand [in3]
    rearCaliperVol = rearCaliper_def*rearPiston_A*rearPiston_num; % rear caliper volume demand [in3]

    frontMC_delta_caliper = frontCaliperVol/frontMC_A; % front master cylinder change in length due to caliper deformation [in]
    rearMC_delta_caliper = rearCaliperVol/rearMC_A; % rear master cylinder change in length due to caliper deformation [in]

    frontMC_delta = frontMC_delta_caliper + biasBar_def + mcMount_def; % front master cylinder change in lengthh due to caliper and bias bar deformation [in]
    rearMC_delta = rearMC_delta_caliper + biasBar_def + mcMount_def; % rear master cylinder change in length due to caliper and bias bar deformation [in]
end