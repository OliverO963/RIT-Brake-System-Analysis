function P = loadVehicleParams(xlsxPath)
%LOADVEHICLEPARAMS Load centralized vehicle/brake parameters from spreadsheet.
%
%   P = loadVehicleParams()      reads VehicleParameters.xlsx from the same
%                                folder as this function.
%   P = loadVehicleParams(path)  reads a specific spreadsheet.
%
%   Returns a struct with one field per row of the "Parameters" sheet, plus
%   the mu-vs-temperature lookup vectors from the "MuTable" sheet, plus a
%   small set of DERIVED convenience fields (see below).
%
%   ============ BASIS CONVENTION ============
%   Every sizing parameter in the spreadsheet is declared strictly
%   PER SINGLE COMPONENT. Specifically:
%       A_pad_front/rear      - ONE pad   (x2 = one corner, x4 = one axle)
%       RotorArea_front/rear  - ONE rotor
%       RotorMass_front/rear  - ONE rotor
%       I_corner              - ONE corner of the vehicle
%       *_piston_count        - ONE caliper, ALL pistons both sides combined
%       *_piston_dia          - ONE piston
%       *_pad_mean_radius     - ONE rotor (friction lever arm)
%   VehicleMass and the aero polynomial are the only whole-vehicle values.
%
%   Consumers must scale UP from these, never assume a coarser basis. The
%   derived fields below exist so the two most error-prone scalings are done
%   in exactly one place:
%
%       P.front_piston_area_per_side / P.rear_piston_area_per_side
%           Total bore area pushing from ONE side of the caliper. For an
%           opposed (fixed) caliper the pistons split evenly across both
%           sides, so only half of them generate clamp force - the other
%           half react against it. The standard disc torque identity
%           T = 2*mu*F_clamp*r_eff already accounts for the two friction
%           faces, so F_clamp must be the SINGLE-SIDE force or the torque
%           comes out 2x high.
%
%       P.A_pad_front_cm2 / P.A_pad_rear_cm2
%           Pad area in cm^2, still per SINGLE pad.

if nargin < 1 || isempty(xlsxPath)
    xlsxPath = fullfile(fileparts(mfilename('fullpath')), 'VehicleParameters.xlsx');
end
if ~isfile(xlsxPath)
    error('loadVehicleParams:FileNotFound', ...
        'Vehicle parameter spreadsheet not found: %s', xlsxPath);
end

% ---- Scalar parameters ----
T = readtable(xlsxPath, 'Sheet', 'Parameters', 'TextType', 'string');
requiredCols = {'Name', 'Value'};
missingCols = requiredCols(~ismember(requiredCols, T.Properties.VariableNames));
if ~isempty(missingCols)
    error('loadVehicleParams:BadSchema', ...
        'Sheet "Parameters" is missing required column(s): %s', strjoin(missingCols, ', '));
end

P = struct();
for i = 1:height(T)
    name = strtrim(char(T.Name(i)));
    if isempty(name)
        continue
    end
    if ~isvarname(name)
        error('loadVehicleParams:BadName', ...
            'Row %d: "%s" is not a valid MATLAB field name.', i, name);
    end
    value = T.Value(i);
    if iscell(value); value = value{1}; end
    if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
        error('loadVehicleParams:BadValue', ...
            'Parameter "%s" (row %d) must be a finite numeric scalar.', name, i);
    end
    if isfield(P, name)
        error('loadVehicleParams:DuplicateName', ...
            'Parameter "%s" is defined more than once.', name);
    end
    P.(name) = value;
end

% ---- Mu vs temperature lookup ----
Tmu = readtable(xlsxPath, 'Sheet', 'MuTable');
if ~all(ismember({'Temp_F', 'Mu'}, Tmu.Properties.VariableNames))
    error('loadVehicleParams:BadSchema', ...
        'Sheet "MuTable" must have columns Temp_F and Mu.');
end
P.mu_temp_table = Tmu.Temp_F(:)';
P.mu_table      = Tmu.Mu(:)';
if numel(P.mu_temp_table) < 2
    error('loadVehicleParams:BadMuTable', 'MuTable needs at least 2 points.');
end

% ---- Presence check: fail loudly on a renamed/deleted row rather than
%      silently propagating an undefined field into the physics ----
expected = {'VehicleMass','I_corner','WheelR','gear_ratio','TambC_fallback', ...
    'RotorMass_front','RotorMass_rear','RotorArea_front','RotorArea_rear', ...
    'front_rotor_dia','rear_rotor_dia','A_pad_front','A_pad_rear','pads_per_corner', ...
    'front_pad_mean_radius','rear_pad_mean_radius', ...
    'front_piston_count','rear_piston_count','front_piston_dia','rear_piston_dia', ...
    'front_caliper_opposed','rear_caliper_opposed','BrakeFrac','CalibrationFactor', ...
    'aero_open_a','aero_open_b','aero_open_c','aero_closed_a','aero_closed_b','aero_closed_c', ...
    'press_front_slope','press_front_offset','press_rear_slope','press_rear_offset', ...
    'temp_adc_slope','temp_adc_offset','velx_threshold','min_pressure', ...
    'current_limit_actual','torque_limit_total_MN','regen_power_limit', ...
    'regen_max_curve_a','regen_max_curve_b','regen_efficiency', ...
    'min_regen_speed','min_energy_threshold'};
missingParams = expected(~isfield(P, expected));
if ~isempty(missingParams)
    error('loadVehicleParams:MissingParams', ...
        'Spreadsheet is missing required parameter(s): %s', strjoin(missingParams, ', '));
end

% ---- Derived convenience fields (see header for why these live here) ----
frontSides = 1 + (P.front_caliper_opposed ~= 0);   % opposed -> 2 sides, floating -> 1
rearSides  = 1 + (P.rear_caliper_opposed  ~= 0);
P.front_piston_area_per_side = (P.front_piston_count / frontSides) * pi * (P.front_piston_dia/2)^2;
P.rear_piston_area_per_side  = (P.rear_piston_count  / rearSides)  * pi * (P.rear_piston_dia/2)^2;

P.A_pad_front_cm2 = P.A_pad_front / 100;   % mm^2 -> cm^2, still ONE pad
P.A_pad_rear_cm2  = P.A_pad_rear  / 100;

P.sourceFile = xlsxPath;
end
