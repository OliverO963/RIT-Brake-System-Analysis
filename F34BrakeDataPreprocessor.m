function manifest = F34BrakeDataPreprocessor(inputFiles, outputFolder)
%F34BRAKEDATAPREPROCESSOR Import, merge, trim, split, and reformat F34 data.
%
%   manifest = F34BrakeDataPreprocessor(inputFiles, outputFolder)
%
%   The function reads one or more selected .txt, .csv, .xlsx, and .xls files,
%   recognizes channels from their headers, merges complementary exports
%   belonging to the same session, removes long stationary periods, splits
%   the session into driving events, removes the inactive intervals, and
%   writes one combined file per session compatible with the existing
%   BrakeCoeffOptimizer.m parser.
%
%   Output layouts:
%       21 columns: BrakeCoeffOptimizer "no_drs" format
%       22 columns: BrakeCoeffOptimizer "drs" format
%
%   The generated text files always contain exactly six header lines because
%   BrakeCoeffOptimizer uses NumHeaderLines = 6.
%
%   Important behavior:
%   - Files with the same normalized base name are merged by time. For
%     example, "Endurance 4-17.xlsx" and "Endurance 4-17(1).txt" are treated
%     as complementary exports from one session.
%   - The script does not invent missing sensor data. Sessions with missing
%     optimizer channels are still exported in the fixed 21/22-column layout,
%     with NaN written in positions corresponding to channels that were not
%     logged. The manifest marks these files as exported_incomplete.
%   - By default, all detected driving segments from a session are concatenated
%     into one output file. Removed idle gaps are compressed out of the output
%     time vector. Set cfg.CombineSegmentsIntoOneFile = false in defaultConfig
%     to restore one output file per detected segment.
%   - Motor RPM may be reconstructed from vehicle speed when RPM channels are
%     absent. Longitudinal acceleration may be reconstructed from speed.
%   - Startup temperature spikes are cleaned with a robust moving-median filter.
%     Data before the first stable temperature window is trimmed, and output time
%     is reset automatically. Files without a stable start can be skipped.
%
%   Examples:
%       manifest = F34BrakeDataPreprocessor;  % opens a file picker
%       manifest = F34BrakeDataPreprocessor("C:\\Data\\Endurance.txt");
%       manifest = F34BrakeDataPreprocessor("C:\\Data\\Endurance.txt", ...
%           "C:\\Data\\Processed");
%       manifest = F34BrakeDataPreprocessor(["C:\\Data\\Endurance.txt", ...
%           "C:\\Data\\Endurance.xlsx"]);
%
%   Tested design target: MATLAB R2021b or newer.

% When no files are supplied, ask the user to choose one or more files.
% Hold Ctrl (Windows) or Command (macOS) to select complementary exports,
% such as the TXT and XLSX files belonging to the same driving session.
if nargin < 1 || isempty(inputFiles)
    fileFilter = { ...
        '*.txt;*.csv;*.xlsx;*.xls', ...
        'F34 data files (*.txt, *.csv, *.xlsx, *.xls)'; ...
        '*.*', 'All files (*.*)'};

    [selectedFiles, selectedPath] = uigetfile(fileFilter, ...
        'Select one or more F34 data files', ...
        'MultiSelect', 'on');

    if isequal(selectedFiles, 0)
        fprintf('File selection canceled. No files were processed.\n');
        manifest = table();
        return;
    end

    if ischar(selectedFiles)
        selectedFiles = {selectedFiles};
    end

    files = string(fullfile(selectedPath, selectedFiles(:)));
    inputFolder = string(selectedPath);
else
    files = string(inputFiles(:));
    if isempty(files)
        fprintf('No input files were supplied.\n');
        manifest = table();
        return;
    end

    [firstFolder, ~, ~] = fileparts(files(1));
    if strlength(firstFolder) == 0
        firstFolder = pwd;
        files = string(fullfile(firstFolder, files));
    end
    inputFolder = string(firstFolder);
end

cfg = defaultConfig();

supportedExtensions = [".txt", ".csv", ".xlsx", ".xls"];
for k = 1:numel(files)
    if ~isfile(files(k))
        error('F34BrakeDataPreprocessor:FileNotFound', ...
            'Input file does not exist: %s', files(k));
    end

    [~, ~, extension] = fileparts(files(k));
    if ~any(strcmpi(string(extension), supportedExtensions))
        error('F34BrakeDataPreprocessor:UnsupportedFile', ...
            'Unsupported input type for %s. Select TXT, CSV, XLSX, or XLS.', ...
            files(k));
    end
end
files = unique(files, 'stable');

% Ask where processed files should be written. No automatic subfolder is
% created when the function is started without an explicit outputFolder.
if nargin < 2 || isempty(outputFolder)
    selectedOutputFolder = uigetdir(char(inputFolder), ...
        'Select the output folder for processed files');

    if isequal(selectedOutputFolder, 0)
        fprintf('Output folder selection canceled. No files were processed.\n');
        manifest = table();
        return;
    end

    outputFolder = string(selectedOutputFolder);
else
    outputFolder = string(outputFolder);
    if ~isfolder(outputFolder)
        mkdir(outputFolder);
    end
end

fprintf('F34 Brake Data Preprocessor\n');
fprintf('Selected files: %d\n', numel(files));
fprintf('Output folder: %s\n\n', outputFolder);

sources = struct([]);
for k = 1:numel(files)
    fprintf('[%d/%d] Reading %s\n', k, numel(files), files(k));
    try
        src = readSourceFile(files(k), cfg);
        if isempty(src.time)
            warning('No usable numeric rows were found in %s. Skipping.', files(k));
            continue;
        end
        if isempty(sources)
            sources = src;
        else
            sources(end+1) = src; %#ok<AGROW>
        end
        fprintf('       %d rows, %d recognized channels, session key "%s"\n', ...
            numel(src.time), src.recognizedCount, src.sessionKey);
        fprintf(['       Source time range: %.3f to %.3f s; temperature samples ', ...
            'FL=%d, FR=%d, RL=%d, RR=%d\n'], ...
            src.time(1), src.time(end), ...
            finiteCount(src.channels.tempFLRaw), ...
            finiteCount(src.channels.tempFRRaw), ...
            finiteCount(src.channels.tempRLRaw), ...
            finiteCount(src.channels.tempRRRaw));
    catch ME
        warning('Could not read %s: %s', files(k), ME.message);
    end
end

if isempty(sources)
    error('F34BrakeDataPreprocessor:NoUsableSources', ...
        'None of the discovered files contained usable time-series data.');
end

sessionKeys = unique(string({sources.sessionKey}), 'stable');
manifestRows = struct([]);

for g = 1:numel(sessionKeys)
    key = sessionKeys(g);
    groupSources = sources(string({sources.sessionKey}) == key);
    fprintf('\nProcessing session "%s" from %d file(s).\n', key, numel(groupSources));

    session = mergeSessionSources(groupSources, cfg);
    session.kind = classifySession(session, cfg);
    [session, readiness] = prepareOptimizerChannels(session, cfg);

    if cfg.RawFormatOnly
        fprintf('  Raw format-only mode: skipping temperature-startup cleaning and segment detection.\n');
    
        startupInfo = struct('stableStartFound', false, 'spikesReplaced', 0, ...
            'trimmedDuration_s', 0, 'startTemperature_C', NaN, ...
            'selectedSensors', strings(0,1), 'rejectedSensors', strings(0,1), ...
            'reason', "skipped_raw_format_only");
        readiness = evaluateOptimizerReadiness(session, cfg);
    
        if cfg.RequireTemperatureData && ~readiness.hasTemperatureData
            fprintf('  Not exported: usable front and rear rotor-temperature data are required.\n');
            manifestRows = appendManifestRow(manifestRows, session, readiness, ...
                0, NaN, NaN, NaN, "", "skipped_missing_temperature", 0);
            continue;
        end
    
        wholeSegment = struct('startIndex', 1, 'endIndex', numel(session.time));
        [outMatrix, outHeaders, formatName] = buildOptimizerMatrix(session, wholeSegment, cfg);
    
        safeKey = sanitizeFileName(key);
        fileName = sprintf('%s_rawformat.txt', safeKey);
        outputPath = fullfile(outputFolder, fileName);
    
        metadata = struct();
        metadata.sessionKey = key;
        metadata.kind = "raw_format_only";
        metadata.startTime = session.time(wholeSegment.startIndex);
        metadata.endTime = session.time(wholeSegment.endIndex);
        metadata.segmentCount = 1;
        metadata.sourceRanges = sprintf('%.3f-%.3f s', metadata.startTime, metadata.endTime);
        metadata.formatName = formatName;
        metadata.sourceFiles = string({groupSources.fileName});
        metadata.estimatedChannels = readiness.estimated;
        metadata.missingChannels = readiness.missing;
        metadata.startupInfo = startupInfo;
    
        writeOptimizerTextFile(outputPath, outMatrix, outHeaders, metadata);
        fprintf('  Wrote %s (%d rows, %d columns) - no trimming or filtering applied.\n', ...
            fileName, size(outMatrix,1), size(outMatrix,2));
    
        exportStatus = "exported_raw_format";
        if ~readiness.ready
            exportStatus = "exported_incomplete_raw_format";
        end
        manifestRows = appendManifestRow(manifestRows, session, readiness, ...
            1, metadata.startTime, metadata.endTime, ...
            metadata.endTime - metadata.startTime, outputPath, exportStatus, 1);
    
        continue;   % skip the normal cleaning/segmenting/export path below
    end

    [session, startupInfo] = cleanTemperatureStartup(session, cfg);
    readiness = evaluateOptimizerReadiness(session, cfg);

    if startupInfo.stableStartFound
        fprintf(['  Temperature startup cleaning: removed %.2f s; ', ...
            'stable start approximately %.1f C (%d sensor(s)).\n'], ...
            startupInfo.trimmedDuration_s, startupInfo.startTemperature_C, ...
            startupInfo.validSensorCount);
        if ~isempty(startupInfo.selectedSensors)
            fprintf('  Temperature pair used for startup: %s.\n', ...
                strjoin(startupInfo.selectedSensors, "/"));
        end
        if ~isempty(startupInfo.rejectedSensors)
            fprintf('  Rejected inconsistent temperature channel(s): %s.\n', ...
                strjoin(startupInfo.rejectedSensors, ", "));
        end
        if startupInfo.spikesReplaced > 0
            fprintf('  Replaced %d isolated temperature spike sample(s).\n', ...
                startupInfo.spikesReplaced);
        end
    elseif cfg.TemperatureStartup.Enable
        fprintf('  Temperature startup cleaning: no stable start found (%s).\n', ...
            startupInfo.reason);
    end

    fprintf('  Usable temperature samples after merge/cleaning: front=%d, rear=%d\n', ...
        readiness.frontTemperatureSamples, readiness.rearTemperatureSamples);

    if cfg.RequireTemperatureData && ~readiness.hasTemperatureData
        fprintf(['  Not exported: usable front and rear rotor-temperature ', ...
            'data are required.\n']);
        manifestRows = appendManifestRow(manifestRows, session, readiness, ...
            0, NaN, NaN, NaN, "", "skipped_missing_temperature", 0);
        continue;
    end

    if cfg.TemperatureStartup.Enable && ...
            cfg.TemperatureStartup.RequireStableStart && ...
            ~startupInfo.stableStartFound
        fprintf('  Not exported: a stable initial rotor-temperature window is required.\n');
        manifestRows = appendManifestRow(manifestRows, session, readiness, ...
            0, NaN, NaN, NaN, "", "skipped_unstable_temperature_start", 0);
        continue;
    end

    segments = detectSegments(session, cfg);

    if isempty(segments)
        fprintf('  No driving segments met the configured thresholds.\n');
        manifestRows = appendManifestRow(manifestRows, session, readiness, ...
            0, NaN, NaN, NaN, "", "no_driving_segment", 0);
        continue;
    end

    fprintf('  Detected %d %s segment(s).\n', numel(segments), session.kind);
    segmentSettings = getSegmentSettings(session.kind, cfg);
    fprintf(['  Retaining up to %.1f s before and up to %.1f s after each ', ...
        'detected stint (reaches toward the neighboring stint or file ', ...
        'boundary when available).\n'], ...
        segmentSettings.PrePadding_s, segmentSettings.PostCooldown_s);

    if ~readiness.ready
        fprintf('  Exporting incomplete data. Missing channels will be NaN: %s\n', ...
            strjoin(readiness.missing, ', '));
    end

    safeKey = sanitizeFileName(key);
    typeLabel = sanitizeFileName(session.kind);

    if cfg.CombineSegmentsIntoOneFile
        [outMatrix, outHeaders, formatName, combinedInfo] = ...
            buildCombinedOptimizerMatrix(session, segments, cfg);

        fileName = sprintf('%s_%s_combined.txt', safeKey, typeLabel);
        outputPath = fullfile(outputFolder, fileName);

        metadata = struct();
        metadata.sessionKey = key;
        metadata.kind = session.kind;
        metadata.startTime = combinedInfo.firstSourceTime;
        metadata.endTime = combinedInfo.lastSourceTime;
        metadata.segmentCount = numel(segments);
        metadata.sourceRanges = combinedInfo.sourceRanges;
        metadata.formatName = formatName;
        metadata.sourceFiles = string({groupSources.fileName});
        metadata.estimatedChannels = readiness.estimated;
        metadata.missingChannels = readiness.missing;
        metadata.startupInfo = startupInfo;

        writeOptimizerTextFile(outputPath, outMatrix, outHeaders, metadata);
        fprintf('  Wrote one combined file: %s (%d rows, %d columns).\n', ...
            fileName, size(outMatrix,1), size(outMatrix,2));

        if readiness.ready
            exportStatus = "exported_combined";
        else
            exportStatus = "exported_incomplete_combined";
        end
        manifestRows = appendManifestRow(manifestRows, session, readiness, ...
            0, combinedInfo.firstSourceTime, combinedInfo.lastSourceTime, ...
            combinedInfo.outputDuration, outputPath, exportStatus, numel(segments));
    else
        for s = 1:numel(segments)
            seg = segments(s);
            [outMatrix, outHeaders, formatName] = buildOptimizerMatrix(session, seg, cfg);

            fileName = sprintf('%s_%s_%03d.txt', safeKey, typeLabel, s);
            outputPath = fullfile(outputFolder, fileName);

            metadata = struct();
            metadata.sessionKey = key;
            metadata.kind = session.kind;
            metadata.startTime = session.time(seg.startIndex);
            metadata.endTime = session.time(seg.endIndex);
            metadata.segmentCount = 1;
            metadata.sourceRanges = sprintf('%.3f-%.3f s', ...
                metadata.startTime, metadata.endTime);
            metadata.formatName = formatName;
            metadata.sourceFiles = string({groupSources.fileName});
            metadata.estimatedChannels = readiness.estimated;
            metadata.missingChannels = readiness.missing;
            metadata.startupInfo = startupInfo;

            writeOptimizerTextFile(outputPath, outMatrix, outHeaders, metadata);
            fprintf('  Wrote %s (%d rows, %d columns).\n', fileName, ...
                size(outMatrix,1), size(outMatrix,2));

            if readiness.ready
                exportStatus = "exported";
            else
                exportStatus = "exported_incomplete";
            end
            manifestRows = appendManifestRow(manifestRows, session, readiness, ...
                s, metadata.startTime, metadata.endTime, ...
                metadata.endTime-metadata.startTime, outputPath, exportStatus, 1);
        end
    end
end

manifest = makeManifestTable(manifestRows);
manifestPath = fullfile(outputFolder, 'PreprocessManifest.csv');
writetable(manifest, manifestPath);

fprintf('\nFinished. Manifest written to:\n  %s\n', manifestPath);
fprintf('Exported files: %d\n', sum(startsWith(manifest.Status, "exported")));
fprintf('Skipped files : %d\n', sum(startsWith(manifest.Status, "skipped")));
end

%% Configuration
function cfg = defaultConfig()
cfg.OutputFolderName = "BrakeCoeffOptimizer_Ready";

% Export all retained driving portions from each session into one file.
% Set this false to return to one file per detected segment.
cfg.CombineSegmentsIntoOneFile = true;

% When true, skip temperature-startup cleaning, spike filtering, and
% segment detection entirely - just reformat each session's full time
% range into the optimizer's column layout and write it out as-is.
% Useful when the input has already been trimmed/curated upstream.
cfg.RawFormatOnly = true;

% F34 vehicle values used only when motor RPM must be estimated from speed.
cfg.GearRatio = 12.97;
cfg.WheelRadius_m = 0.213;

% Sensor validity limits.
cfg.MaxValidSpeed_mps = 40;
cfg.MaxValidAbsAccel_mps2 = 50;
cfg.MaxInterpolationGap_s = 1.0;

% Automatically reconcile complementary exports that describe the same
% recording but use different time origins, such as one file beginning at
% 0 s and another beginning at 3500 s. Absolute timestamps are retained when
% the files already overlap. A shift is attempted only when their durations
% are similar and relative-time alignment substantially improves overlap.
cfg.AutoAlignTimeOrigins = true;
cfg.MinimumAbsoluteTimeOverlapFraction = 0.25;
cfg.MaximumDurationMismatchFraction = 0.15;

% General movement/activity detection.
cfg.MovingSpeedThreshold_mps = 1.0;
cfg.AccelActivityThreshold_mps2 = 0.35;
cfg.ActivityWindow_s = 1.0;
cfg.ActivityFractionThreshold = 0.15;
cfg.BrakeActivityDelta_ADC = 20;

% Segment settings by session type.
%
% PrePadding_s / PostCooldown_s are now MAXIMUM reach caps, not fixed
% amounts. Each detected stint reaches back toward the end of the
% previous stint (or the start of the file, for the first stint) and
% forward toward the start of the next stint (or the end of the file,
% for the last stint), up to these caps. This means:
%   - A stint normally gets its full available lead-in (cold start) and
%     cooldown data, not just a token few seconds, fixing exports whose
%     t=0 temperature was already several hundred F.
%   - The cooldown tail is only truncated short of the next stint when
%     the idle gap between stints genuinely exceeds the cap (e.g. a
%     lunch break), which is "where sensible" -- not by default.
%   - Two stints separated by a short gap end up with overlapping
%     pre/post windows and are automatically stitched into one
%     continuous file by mergeOverlappingSegments, avoiding an
%     artificial temperature discontinuity at the cut.
cfg.Endurance.MergeIdleGap_s = 30;
cfg.Endurance.MinimumDuration_s = 30;
cfg.Endurance.PrePadding_s = 30;      % max data retained before driving begins
cfg.Endurance.PostCooldown_s = 300;   % max cooldown data retained after driving ends

cfg.Autocross.MergeIdleGap_s = 8;
cfg.Autocross.MinimumDuration_s = 10;
cfg.Autocross.PrePadding_s = 15;
cfg.Autocross.PostCooldown_s = 120;

cfg.Acceleration.MergeIdleGap_s = 4;
cfg.Acceleration.MinimumDuration_s = 2.5;
cfg.Acceleration.PrePadding_s = 5;
cfg.Acceleration.PostCooldown_s = 45;

cfg.Braking.MergeIdleGap_s = 1.0;
cfg.Braking.MinimumDuration_s = 0.5;
cfg.Braking.PrePadding_s = 2.0;
cfg.Braking.PostCooldown_s = 30;
cfg.Braking.MinimumDecel_mps2 = -0.5;
cfg.Braking.MinimumSpeed_mps = 2.0;

cfg.Moving.MergeIdleGap_s = 10;
cfg.Moving.MinimumDuration_s = 3;
cfg.Moving.PrePadding_s = 10;
cfg.Moving.PostCooldown_s = 120;

% Noise-robustness settings used during segment detection. A driving
% stint should show sustained activity, not a couple of isolated
% seconds of it. Without these, a cluster of brief sensor glitches
% (ADC noise, EMI, etc.) spaced less than MergeIdleGap_s apart gets
% bridged into one span that is long enough to pass MinimumDuration_s
% even though almost none of it was real activity.
cfg.NoiseRejection.MinRawBlip_s = 0.3;        % drop isolated blips shorter than this before bridging
cfg.NoiseRejection.MinActivityDensity = 0.20; % bridged run must be >=20% real "active" samples, else rejected

% Require usable front- and rear-rotor temperature data before export.
% BrakeCoeffOptimizer fits measured temperatures, so a file without both
% axle temperatures cannot produce a meaningful coefficient result.
cfg.RequireTemperatureData = true;
cfg.MinimumTemperatureSamples = 10;

% Startup temperature cleaning. Analog temperature sensors can experience
% EMI spikes while the HV system is energizing. The preprocessor replaces
% isolated spikes, locates the first persistent stable temperature window,
% trims earlier rows, and resets the session time without forcing ambient.
cfg.TemperatureStartup.Enable = true;
cfg.TemperatureStartup.RequireStableStart = true;
cfg.TemperatureStartup.StabilityWindow_s = 3.0;
cfg.TemperatureStartup.ConfirmationWindow_s = 3.0;
cfg.TemperatureStartup.MaximumSearchTime_s = 25.0;
cfg.TemperatureStartup.MaxWindowRange_C = 4.0;
cfg.TemperatureStartup.MaxMedianRate_C_per_s = 1.5;
cfg.TemperatureStartup.MaxConfirmationShift_C = 5.0;
cfg.TemperatureStartup.MaxSensorSpread_C = 30.0;
cfg.TemperatureStartup.MaximumSpeed_mps = 1.0;
cfg.TemperatureStartup.MinimumValidSensors = 2;

% The optimizer reads the right-front and right-rear temperature positions.
% Evaluate the right-side pair first, but allow the left-side pair as a
% fallback when the right-side sensors are unavailable.
cfg.TemperatureStartup.PreferredSide = "right";
cfg.TemperatureStartup.AllowOppositeSideFallback = true;
cfg.TemperatureStartup.AllowMixedSideFallback = false;

% If a non-selected corner disagrees strongly with the coherent selected
% pair at startup, treat that corner as an invalid/noisy channel and export
% it as NaN rather than allowing it to block the whole session.
cfg.TemperatureStartup.RejectDisagreeingChannels = true;
cfg.TemperatureStartup.MaxInitialChannelDisagreement_C = 60;

cfg.TemperatureStartup.MinimumPlausible_C = -30;
cfg.TemperatureStartup.MaximumPlausible_C = 800;

% Keep left/right temperature channels independent. Missing channels are not
% copied across an axle unless this is explicitly enabled.
cfg.CopyMissingTemperatureAcrossAxle = false;

% Robust isolated-spike filter. Applied across the ENTIRE file, not just
% the startup window: a rotor has enough thermal mass that a genuine
% temperature change can't spike and fall back to ambient within a
% couple of seconds, so an isolated spike-then-return pattern anywhere
% in the file is sensor/EMI noise, not real physics. The filter compares
% each sample against a local robust median (movmedian) and only
% replaces points that deviate more than a MAD-based threshold, so a
% genuine multi-sample heating RAMP (which stays elevated, not an
% isolated point) is left untouched.
cfg.TemperatureStartup.SpikeFilter.Enable = true;
cfg.TemperatureStartup.SpikeFilter.ApplyToEntireFile = true;
cfg.TemperatureStartup.SpikeFilter.Window_s = 1.0;
cfg.TemperatureStartup.SpikeFilter.MinimumDeviation_C = 12.0;
cfg.TemperatureStartup.SpikeFilter.MadFactor = 6.0;

% Safety policy for incomplete data.
% The supplied optimizer uses motor torque to calculate regen energy. Keeping
% this false prevents the script from silently assuming zero regen.
cfg.AllowZeroMotorTorqueFallback = false;
end

%% File discovery and import
function files = findInputFiles(inputFolder, outputFolder)
extensions = ["*.txt", "*.csv", "*.xlsx", "*.xls"];
files = strings(0,1);
for e = 1:numel(extensions)
    listing = dir(fullfile(inputFolder, '**', extensions(e)));
    for k = 1:numel(listing)
        fullName = string(fullfile(listing(k).folder, listing(k).name));
        if startsWith(fullName, string(outputFolder), 'IgnoreCase', true)
            continue;
        end
        files(end+1,1) = fullName; %#ok<AGROW>
    end
end
files = unique(files, 'stable');
end

function src = readSourceFile(filePath, cfg)
[~, fileName, extension] = fileparts(filePath);
extension = lower(string(extension));

switch extension
    case {".txt", ".csv"}
        [headers, data] = readDelimitedFile(filePath);
    case {".xlsx", ".xls"}
        [headers, data] = readSpreadsheetFile(filePath);
    otherwise
        error('Unsupported file type: %s', extension);
end

if isempty(data)
    src = emptySource(filePath, fileName);
    return;
end

[headers, data] = trimUnusedColumns(headers, data);
channels = recognizeChannels(headers, data);

if isempty(channels.time)
    error('No time column was recognized in %s.', filePath);
end

% Sort by time and remove duplicate timestamps.
time = double(channels.time(:));
validTime = isfinite(time);
time = time(validTime);
data = data(validTime,:);
[time, order] = sort(time);
data = data(order,:);
[time, uniqueRows] = unique(time, 'stable');
data = data(uniqueRows,:);
channels = recognizeChannels(headers, data);
channels.time = time;

fields = canonicalFieldNames();
recognizedCount = 0;
for k = 1:numel(fields)
    if ~isempty(channels.(fields(k)))
        recognizedCount = recognizedCount + 1;
    end
end

src = struct();
src.filePath = string(filePath);
src.fileName = string(fileName) + extension;
src.sessionKey = normalizeSessionKey(fileName);
src.headers = headers;
src.data = data;
src.time = time;
src.channels = channels;
src.recognizedCount = recognizedCount;
src.qualityScore = recognizedCount + log10(max(numel(time),1));
src.cfg = cfg;
end

function src = emptySource(filePath, fileName)
src = struct('filePath', string(filePath), 'fileName', string(fileName), ...
    'sessionKey', normalizeSessionKey(fileName), 'headers', strings(1,0), ...
    'data', [], 'time', [], 'channels', struct(), 'recognizedCount', 0, ...
    'qualityScore', 0, 'cfg', struct());
end

function [headers, data] = readDelimitedFile(filePath)
lines = readlines(filePath, 'EmptyLineRule', 'read');
if isempty(lines)
    headers = strings(1,0);
    data = [];
    return;
end

headerRow = findTextHeaderRow(lines);
headerLine = lines(headerRow);
delimiter = detectDelimiter(headerLine);
headers = strtrim(split(headerLine, delimiter)).';

try
    data = readmatrix(filePath, 'FileType', 'text', ...
        'Delimiter', char(delimiter), 'NumHeaderLines', headerRow, ...
        'TreatAsMissing', {'', 'NaN', 'nan', 'N/A'});
catch
    % Fallback for unusual text exports.
    rawLines = lines(headerRow+1:end);
    data = nan(numel(rawLines), numel(headers));
    for r = 1:numel(rawLines)
        pieces = split(rawLines(r), delimiter);
        n = min(numel(pieces), numel(headers));
        data(r,1:n) = str2double(pieces(1:n));
    end
end

if isempty(data)
    return;
end
if isvector(data)
    data = data(:);
end
end

function headerRow = findTextHeaderRow(lines)
headerRow = [];
limit = min(numel(lines), 100);
for r = 1:limit
    line = strtrim(lines(r));
    if startsWith(line, "#") || strlength(line) == 0
        continue;
    end
    normalized = normalizeHeader(line);
    if contains(normalized, "xtime") || startsWith(normalized, "time")
        headerRow = r;
        break;
    end
end
if isempty(headerRow)
    error('Could not locate a time-series header row in the text file.');
end
end

function delimiter = detectDelimiter(headerLine)
candidates = [sprintf('\t'), ',', ';'];
counts = zeros(size(candidates));
for k = 1:numel(candidates)
    counts(k) = count(headerLine, candidates(k));
end
[~, idx] = max(counts);
delimiter = string(candidates(idx));
end

function [headers, data] = readSpreadsheetFile(filePath)
sheets = sheetnames(filePath);
if isempty(sheets)
    error('The workbook contains no readable worksheets.');
end

% Use the first worksheet containing a recognizable time header.
allCells = [];
headerRow = [];
for s = 1:numel(sheets)
    try
        candidate = readcell(filePath, 'Sheet', sheets(s));
    catch
        continue;
    end
    candidateHeader = findCellHeaderRow(candidate);
    if ~isempty(candidateHeader)
        allCells = candidate;
        headerRow = candidateHeader;
        break;
    end
end

if isempty(headerRow)
    error('Could not locate a time-series header row in the workbook.');
end

headers = strtrim(string(allCells(headerRow,:)));
dataCells = allCells(headerRow+1:end,:);
data = str2double(string(dataCells));

% Some spreadsheet readers return numeric cells as missing strings. Patch
% those cells directly from the original cell array.
for c = 1:size(dataCells,2)
    numericMask = cellfun(@isnumeric, dataCells(:,c));
    if any(numericMask)
        values = dataCells(numericMask,c);
        data(numericMask,c) = cellfun(@double, values);
    end
end
end

function headerRow = findCellHeaderRow(cells)
headerRow = [];
if isempty(cells)
    return;
end
limit = min(size(cells,1), 100);
for r = 1:limit
    rowText = string(cells(r,:));
    normalized = normalizeHeader(rowText);
    timeAliases = ["xtimes", "xtime", "times", "time"];
    if any(ismember(normalized, timeAliases)) && sum(strlength(strtrim(rowText)) > 0) >= 2
        headerRow = r;
        return;
    end
end
end

function [headers, data] = trimUnusedColumns(headers, data)
headers = string(headers(:).');
if isempty(data)
    return;
end
n = min(numel(headers), size(data,2));
headers = headers(1:n);
data = data(:,1:n);

used = strlength(strtrim(headers)) > 0 | any(isfinite(data),1);
lastUsed = find(used, 1, 'last');
if isempty(lastUsed)
    headers = strings(1,0);
    data = [];
    return;
end
headers = headers(1:lastUsed);
data = data(:,1:lastUsed);

for c = 1:numel(headers)
    if ismissing(headers(c)) || strlength(strtrim(headers(c))) == 0
        headers(c) = "UnnamedColumn" + c;
    end
end

validRows = any(isfinite(data),2);
data = data(validRows,:);
end

%% Header recognition
function channels = recognizeChannels(headers, data)
normalized = normalizeHeader(headers);
channels = emptyChannels();

channels.time = getChannel(data, findHeader(normalized, ...
    ["xtimes", "xtime", "times", "time"]));

% Temperature channels are converted to the raw ADC scale expected by
% BrakeCoeffOptimizer. Explicit TEMP/temperature names are preferred. On F34,
% the same four analog temperature inputs may instead be exported under the
% legacy names 1_SSDB_suspension_FL/FR/RL/RR_raw[none]; those names are used as
% corner-specific raw-temperature fallbacks. Celsius/Fahrenheit columns are
% converted back using T_C = 0.246*(ADC - 406).
channels.tempFLRaw = recognizeTemperatureAsRaw(normalized, data, ...
    ["flrotortempraw", "rotortempflraw", "flbraketempraw", "braketempflraw", ...
     "flrotortemperatureraw", "ssdbsuspensionflraw", "suspensionflraw"], ...
    ["flrotortempc", "rotortempflc", "flbraketempc", "flrotortemperaturec"], ...
    ["flrotortempf", "rotortempflf", "flbraketempf", "flrotortemperaturef"], "fl");
channels.tempFRRaw = recognizeTemperatureAsRaw(normalized, data, ...
    ["frrotortempraw", "rotortempfrraw", "frbraketempraw", "braketempfrraw", ...
     "frrotortemperatureraw", "ssdbsuspensionfrraw", "suspensionfrraw"], ...
    ["frrotortempc", "rotortempfrc", "frbraketempc", "frrotortemperaturec"], ...
    ["frrotortempf", "rotortempfrf", "frbraketempf", "frrotortemperaturef"], "fr");
channels.tempRLRaw = recognizeTemperatureAsRaw(normalized, data, ...
    ["rlrotortempraw", "rotortemprlraw", "rlbraketempraw", "braketemprlraw", ...
     "rlrotortemperatureraw", "ssdbsuspensionrlraw", "suspensionrlraw"], ...
    ["rlrotortempc", "rotortemprlc", "rlbraketempc", "rlrotortemperaturec"], ...
    ["rlrotortempf", "rotortemprlf", "rlbraketempf", "rlrotortemperaturef"], "rl");
channels.tempRRRaw = recognizeTemperatureAsRaw(normalized, data, ...
    ["rrrotortempraw", "rotortemprrraw", "rrbraketempraw", "braketemprrraw", ...
     "rrrotortemperatureraw", "ssdbsuspensionrrraw", "suspensionrrraw"], ...
    ["rrrotortempc", "rotortemprrc", "rrbraketempc", "rrrotortemperaturec"], ...
    ["rrrotortempf", "rotortemprrf", "rrbraketempf", "rrrotortemperaturef"], "rr");
channels.tempFrontRaw = recognizeTemperatureAsRaw(normalized, data, ...
    ["frontrotortempraw", "frontbraketempraw", "rotortempfrontraw", ...
     "frontrotortemperatureraw"], ...
    ["frontrotortempc", "frontbraketempc", "rotortempfrontc", ...
     "frontrotortemperaturec"], ...
    ["frontrotortempf", "frontbraketempf", "rotortempfrontf", ...
     "frontrotortemperaturef"], "front");
channels.tempRearRaw = recognizeTemperatureAsRaw(normalized, data, ...
    ["rearrotortempraw", "rearbraketempraw", "rotortemprearraw", ...
     "rearrotortemperatureraw"], ...
    ["rearrotortempc", "rearbraketempc", "rotortemprearc", ...
     "rearrotortemperaturec"], ...
    ["rearrotortempf", "rearbraketempf", "rotortemprearf", ...
     "rearrotortemperaturef"], "rear");

channels.frontPressureRaw = getChannel(data, findHeader(normalized, ...
    ["ssdbbrakepressurefrontraw", "brakepressurefrontraw", ...
     "frontbrakepressureraw", "frontbpraw", "frontbrakepressureadc"], ...
    ["psi", "bar", "kpa", "mpa"]));
channels.rearPressureRaw = getChannel(data, findHeader(normalized, ...
    ["pedalinputsrawbrakesrearadc", "brakesrearadc", ...
     "rearbrakepressureraw", "rearbpraw", "rearbrakepressureadc"], ...
    ["psi", "bar", "kpa", "mpa"]));

channels.steerPct = getChannel(data, findHeader(normalized, ...
    ["inputssteerpct", "steerpct", "steeringpct", "steeringangle"]));
channels.accelPedalA = getChannel(data, findHeader(normalized, ...
    ["accelpositiona", "acceleratorpositiona", "appedala"]));
channels.accelPedalB = getChannel(data, findHeader(normalized, ...
    ["accelpositionb", "acceleratorpositionb", "appedalb"]));

channels.torqueFL = getChannel(data, findHeader(normalized, ...
    ["flfeedbacktorque", "flmotortorque", "torquefl"]));
channels.rpmFL = getChannel(data, findHeader(normalized, ...
    ["flfeedbackvelocityrpm", "flfeedbackvelocity", "flmotorvelocity", "rpmfl"]));
channels.torqueFR = getChannel(data, findHeader(normalized, ...
    ["frfeedbacktorque", "frmotortorque", "torquefr"]));
channels.rpmFR = getChannel(data, findHeader(normalized, ...
    ["frfeedbackvelocityrpm", "frfeedbackvelocity", "frmotorvelocity", "rpmfr"]));
channels.torqueRL = getChannel(data, findHeader(normalized, ...
    ["rlfeedbacktorque", "rlmotortorque", "torquerl"]));
channels.rpmRL = getChannel(data, findHeader(normalized, ...
    ["rlfeedbackvelocityrpm", "rlfeedbackvelocity", "rlmotorvelocity", "rpmrl"]));
channels.torqueRR = getChannel(data, findHeader(normalized, ...
    ["rrfeedbacktorque", "rrmotortorque", "torquerr"]));
channels.rpmRR = getChannel(data, findHeader(normalized, ...
    ["rrfeedbackvelocityrpm", "rrfeedbackvelocity", "rrmotorvelocity", "rpmrr"]));

channels.speed = getChannel(data, findHeader(normalized, ...
    ["vectornavvelbodyxms", "vectornavvelbodyx", "velocityms", ...
     "vehiclespeedms", "vehiclespeed", "velx"], ...
    ["feedbackvelocity", "motorvelocity", "rpm"]));
channels.accelX = getChannel(data, findHeader(normalized, ...
    ["accelxcalibrated", "accelms2", "longitudinalaccel", ...
     "vectornavaccelbodyx", "accelerationx"], ...
    ["position", "pedal"]));
channels.drsState = getChannel(data, findHeader(normalized, ...
    ["drsstate", "drsposition", "drscommand"]));
end

function channels = emptyChannels()
fields = canonicalFieldNames();
for k = 1:numel(fields)
    channels.(fields(k)) = [];
end
end

function fields = canonicalFieldNames()
fields = ["time", "tempFLRaw", "tempFRRaw", "tempRLRaw", "tempRRRaw", ...
    "tempFrontRaw", "tempRearRaw", "frontPressureRaw", "rearPressureRaw", ...
    "steerPct", "accelPedalA", "accelPedalB", ...
    "torqueFL", "rpmFL", "torqueFR", "rpmFR", ...
    "torqueRL", "rpmRL", "torqueRR", "rpmRR", ...
    "speed", "accelX", "drsState"];
end

function idx = findHeader(normalizedHeaders, includePatterns, excludePatterns)
if nargin < 3
    excludePatterns = strings(0,1);
end
idx = [];
for p = 1:numel(includePatterns)
    for h = 1:numel(normalizedHeaders)
        candidate = normalizedHeaders(h);
        if contains(candidate, includePatterns(p))
            excluded = false;
            for e = 1:numel(excludePatterns)
                if contains(candidate, excludePatterns(e))
                    excluded = true;
                    break;
                end
            end
            if ~excluded
                idx = h;
                return;
            end
        end
    end
end
end

function value = getChannel(data, idx)
if isempty(idx) || idx > size(data,2)
    value = [];
else
    value = double(data(:,idx));
end
end

function raw = recognizeTemperatureAsRaw(normalizedHeaders, data, ...
        rawPatterns, celsiusPatterns, fahrenheitPatterns, positionCode)
% Return temperature on the ADC scale used by BrakeCoeffOptimizer.
% Explicit aliases are checked first. A token-based fallback then recognizes
% longer Bosch/DARAB names such as FrontRightBrakeRotorTemperatureRAW.
rawIdx = findHeader(normalizedHeaders, rawPatterns);
if isempty(rawIdx)
    rawIdx = findTemperatureHeader(normalizedHeaders, positionCode, "raw");
end
raw = getChannel(data, rawIdx);
if ~isempty(raw)
    return;
end

cIdx = findHeader(normalizedHeaders, celsiusPatterns);
if isempty(cIdx)
    cIdx = findTemperatureHeader(normalizedHeaders, positionCode, "celsius");
end
tempC = getChannel(data, cIdx);
if ~isempty(tempC)
    raw = tempC ./ 0.246 + 406;
    return;
end

fIdx = findHeader(normalizedHeaders, fahrenheitPatterns);
if isempty(fIdx)
    fIdx = findTemperatureHeader(normalizedHeaders, positionCode, "fahrenheit");
end
tempF = getChannel(data, fIdx);
if ~isempty(tempF)
    tempC = (tempF - 32) .* (5/9);
    raw = tempC ./ 0.246 + 406;
    return;
end

raw = [];
end

function idx = findTemperatureHeader(normalizedHeaders, positionCode, unitKind)
idx = [];
for h = 1:numel(normalizedHeaders)
    candidate = normalizedHeaders(h);
    if ~(contains(candidate, "temp") || contains(candidate, "temperature"))
        continue;
    end
    if contains(candidate, "pressure") || contains(candidate, "ambient") || ...
            contains(candidate, "motor") || contains(candidate, "inverter")
        continue;
    end
    if ~temperaturePositionMatches(candidate, positionCode)
        continue;
    end
    if ~temperatureUnitMatches(candidate, unitKind)
        continue;
    end
    idx = h;
    return;
end
end

function tf = temperaturePositionMatches(candidate, positionCode)
switch string(positionCode)
    case "fl"
        aliases = ["frontleft", "leftfront", "flrotor", "rotorfl", ...
            "flbrake", "brakefl", "fltemp", "tempfl", ...
            "fltemperature", "temperaturefl"];
    case "fr"
        aliases = ["frontright", "rightfront", "frrotor", "rotorfr", ...
            "frbrake", "brakefr", "frtemp", "tempfr", ...
            "frtemperature", "temperaturefr"];
    case "rl"
        aliases = ["rearleft", "leftrear", "rlrotor", "rotorrl", ...
            "rlbrake", "brakerl", "rltemp", "temprl", ...
            "rltemperature", "temperaturerl"];
    case "rr"
        aliases = ["rearright", "rightrear", "rrrotor", "rotorrr", ...
            "rrbrake", "brakerr", "rrtemp", "temprr", ...
            "rrtemperature", "temperaturerr"];
    case "front"
        aliases = ["frontrotor", "rotorfront", "frontbrake", "brakefront", ...
            "fronttemp", "tempfront", "fronttemperature", "temperaturefront"];
    case "rear"
        aliases = ["rearrotor", "rotorrear", "rearbrake", "brakerear", ...
            "reartemp", "temprear", "reartemperature", "temperaturerear"];
    otherwise
        aliases = strings(0,1);
end
tf = any(contains(candidate, aliases));
end

function tf = temperatureUnitMatches(candidate, unitKind)
switch string(unitKind)
    case "raw"
        tf = contains(candidate, "raw") || contains(candidate, "adc");
    case "celsius"
        tf = contains(candidate, "celsius") || contains(candidate, "degc") || ...
            endsWith(candidate, "tempc") || endsWith(candidate, "temperaturec");
    case "fahrenheit"
        tf = contains(candidate, "fahrenheit") || contains(candidate, "degf") || ...
            endsWith(candidate, "tempf") || endsWith(candidate, "temperaturef");
    otherwise
        tf = false;
end
end

function count = finiteCount(value)
if isempty(value)
    count = 0;
else
    count = sum(isfinite(value));
end
end

function normalized = normalizeHeader(text)
normalized = lower(string(text));
normalized = regexprep(normalized, '[^a-z0-9]+', '');
end

%% Session grouping and merging
function key = normalizeSessionKey(fileName)
key = lower(string(fileName));
key = regexprep(key, '\(\s*\d+\s*\)$', '');
key = regexprep(key, '\bcopy\b', '');
key = regexprep(key, '[^a-z0-9]+', '');
if strlength(key) == 0
    key = "session";
end
end

function session = mergeSessionSources(groupSources, cfg)
[~, refIndex] = max([groupSources.qualityScore]);
reference = groupSources(refIndex);
referenceTime = reference.time(:);

% Resolve each source's time origin once, before aligning individual channels.
% This preserves absolute time when files already overlap and shifts only
% complementary exports whose relative timelines clearly match better.
alignedSourceTimes = cell(numel(groupSources),1);
for s = 1:numel(groupSources)
    [alignedSourceTimes{s}, appliedShift] = resolveTimeOrigin( ...
        groupSources(s).time, referenceTime, cfg);
    if abs(appliedShift) > 1e-9
        fprintf('  Time origin adjusted for %s by %+.3f s.\n', ...
            groupSources(s).fileName, appliedShift);
    end
end

session = struct();
session.key = reference.sessionKey;
session.time = referenceTime;
session.sourceFiles = string({groupSources.fileName});
session.channels = emptyChannels();
session.channelSources = struct();
session.estimatedChannels = strings(0,1);
session.cfg = cfg;

fields = canonicalFieldNames();
fields(fields == "time") = [];

for f = 1:numel(fields)
    field = fields(f);
    bestValue = [];
    bestSource = "";
    bestScore = -Inf;

    for s = 1:numel(groupSources)
        candidate = groupSources(s).channels.(field);
        if isempty(candidate)
            continue;
        end
        aligned = alignChannel(alignedSourceTimes{s}, candidate, referenceTime, field);
        finiteFraction = mean(isfinite(aligned));
        score = finiteFraction * 100 + groupSources(s).qualityScore;
        if score > bestScore
            bestValue = aligned;
            bestSource = groupSources(s).fileName;
            bestScore = score;
        end
    end

    session.channels.(field) = bestValue;
    session.channelSources.(field) = bestSource;
end
session.channels.time = referenceTime;
end

function [mappedTime, appliedShift] = resolveTimeOrigin(sourceTime, referenceTime, cfg)
sourceTime = double(sourceTime(:));
referenceTime = double(referenceTime(:));
mappedTime = sourceTime;
appliedShift = 0;

if ~cfg.AutoAlignTimeOrigins || isempty(sourceTime) || isempty(referenceTime)
    return;
end

absoluteOverlap = timeOverlapFraction(sourceTime, referenceTime);
if absoluteOverlap >= cfg.MinimumAbsoluteTimeOverlapFraction
    return;
end

sourceSpan = max(sourceTime) - min(sourceTime);
referenceSpan = max(referenceTime) - min(referenceTime);
if sourceSpan <= 0 || referenceSpan <= 0
    return;
end

durationMismatch = abs(sourceSpan-referenceSpan) / max(sourceSpan, referenceSpan);
if durationMismatch > cfg.MaximumDurationMismatchFraction
    return;
end

candidateShift = referenceTime(1) - sourceTime(1);
candidateTime = sourceTime + candidateShift;
relativeOverlap = timeOverlapFraction(candidateTime, referenceTime);

if relativeOverlap > absoluteOverlap + 0.25 && ...
        relativeOverlap >= cfg.MinimumAbsoluteTimeOverlapFraction
    mappedTime = candidateTime;
    appliedShift = candidateShift;
end
end

function fraction = timeOverlapFraction(timeA, timeB)
a0 = min(timeA); a1 = max(timeA);
b0 = min(timeB); b1 = max(timeB);
overlap = max(0, min(a1,b1) - max(a0,b0));
shorterSpan = min(a1-a0, b1-b0);
if shorterSpan <= 0
    fraction = 0;
else
    fraction = overlap / shorterSpan;
end
end

function aligned = alignChannel(sourceTime, sourceValue, referenceTime, field)
sourceTime = double(sourceTime(:));
sourceValue = double(sourceValue(:));
valid = isfinite(sourceTime) & isfinite(sourceValue);
sourceTime = sourceTime(valid);
sourceValue = sourceValue(valid);

if isempty(sourceTime)
    aligned = nan(size(referenceTime));
    return;
end

[sourceTime, ia] = unique(sourceTime, 'stable');
sourceValue = sourceValue(ia);

if numel(sourceTime) == numel(referenceTime) && ...
        max(abs(sourceTime-referenceTime), [], 'omitnan') < 1e-8
    aligned = sourceValue;
    return;
end

if numel(sourceTime) == 1
    aligned = repmat(sourceValue, size(referenceTime));
elseif field == "drsState"
    aligned = interp1(sourceTime, sourceValue, referenceTime, 'nearest', NaN);
else
    aligned = interp1(sourceTime, sourceValue, referenceTime, 'linear', NaN);
end
end

function kind = classifySession(session, ~)
text = lower(strjoin([session.key; session.sourceFiles(:)], " "));
if contains(text, "endurance")
    kind = "endurance_stint";
elseif contains(text, "autox") || contains(text, "autocross") || contains(text, "skidpad")
    kind = "autocross_run";
elseif contains(text, "accel") || contains(text, "acceleration")
    kind = "acceleration_run";
elseif contains(text, "brake") || contains(text, "braking") || contains(text, "stopping")
    kind = "braking_event";
else
    kind = "moving_event";
end
end

%% Channel preparation and validation
function [session, readiness] = prepareOptimizerChannels(session, cfg)
c = session.channels;
t = session.time(:);
estimated = strings(0,1);

% Reject logger overflow values before interpolation or event detection.
c.speed = invalidateOutside(c.speed, -2, cfg.MaxValidSpeed_mps);
c.accelX = invalidateOutside(c.accelX, -cfg.MaxValidAbsAccel_mps2, cfg.MaxValidAbsAccel_mps2);

% Expand axle temperature measurements to the four positions expected by
% BrakeCoeffOptimizer. The current optimizer uses FR and RR temperatures,
% but all four columns are kept valid for format compatibility.
if isempty(c.tempFLRaw) && ~isempty(c.tempFrontRaw)
    c.tempFLRaw = c.tempFrontRaw;
    estimated(end+1) = "tempFLRaw_from_front_axle"; %#ok<AGROW>
end
if isempty(c.tempFRRaw) && ~isempty(c.tempFrontRaw)
    c.tempFRRaw = c.tempFrontRaw;
    estimated(end+1) = "tempFRRaw_from_front_axle"; %#ok<AGROW>
end
if isempty(c.tempRLRaw) && ~isempty(c.tempRearRaw)
    c.tempRLRaw = c.tempRearRaw;
    estimated(end+1) = "tempRLRaw_from_rear_axle"; %#ok<AGROW>
end
if isempty(c.tempRRRaw) && ~isempty(c.tempRearRaw)
    c.tempRRRaw = c.tempRearRaw;
    estimated(end+1) = "tempRRRaw_from_rear_axle"; %#ok<AGROW>
end

% Keep left and right channels independent by default. Copying a valid
% corner across an axle can hide a failed sensor and makes diagnostics harder.
if cfg.CopyMissingTemperatureAcrossAxle
    [c.tempFLRaw, c.tempFRRaw, estimated] = completeAxlePair( ...
        c.tempFLRaw, c.tempFRRaw, "tempFLRaw", "tempFRRaw", estimated);
    [c.tempRLRaw, c.tempRRRaw, estimated] = completeAxlePair( ...
        c.tempRLRaw, c.tempRRRaw, "tempRLRaw", "tempRRRaw", estimated);
end

% Reconstruct speed from motor RPM when possible.
rpmFields = ["rpmFL", "rpmFR", "rpmRL", "rpmRR"];
if isempty(c.speed)
    rpmMatrix = collectAvailableChannels(c, rpmFields, numel(t));
    if ~isempty(rpmMatrix)
        wheelRPM = median(abs(rpmMatrix), 2, 'omitnan') / cfg.GearRatio;
        c.speed = wheelRPM * (2*pi/60) * cfg.WheelRadius_m;
        estimated(end+1) = "speed_from_motor_rpm"; %#ok<AGROW>
    end
end

% Reconstruct motor RPM from speed when needed. This preserves wheel
% rotational energy calculations but does not invent regen torque.
if ~isempty(c.speed)
    estimatedMotorRPM = c.speed / cfg.WheelRadius_m * (60/(2*pi)) * cfg.GearRatio;
    for k = 1:numel(rpmFields)
        field = rpmFields(k);
        if isempty(c.(field))
            c.(field) = estimatedMotorRPM;
            estimated(end+1) = field + "_from_speed"; %#ok<AGROW>
        end
    end
end

% Derive acceleration from speed if it was not logged.
if isempty(c.accelX) && ~isempty(c.speed)
    speedForDerivative = fillAllMissing(c.speed, t);
    c.accelX = gradient(speedForDerivative, t);
    estimated(end+1) = "accelX_from_speed"; %#ok<AGROW>
end

% Optional zero-regen assumption. Disabled by default because the optimizer
% explicitly uses motor torque to subtract regenerated energy.
torqueFields = ["torqueFL", "torqueFR", "torqueRL", "torqueRR"];
if cfg.AllowZeroMotorTorqueFallback
    for k = 1:numel(torqueFields)
        field = torqueFields(k);
        if isempty(c.(field))
            c.(field) = zeros(size(t));
            estimated(end+1) = field + "_assumed_zero"; %#ok<AGROW>
        end
    end
end

% Nonessential layout columns may safely be zero-filled.
optionalZeroFields = ["steerPct", "accelPedalA", "accelPedalB"];
for k = 1:numel(optionalZeroFields)
    field = optionalZeroFields(k);
    if isempty(c.(field))
        c.(field) = zeros(size(t));
    end
end

required = ["tempFLRaw", "tempFRRaw", "tempRLRaw", "tempRRRaw", ...
    "frontPressureRaw", "rearPressureRaw", ...
    "torqueFL", "rpmFL", "torqueFR", "rpmFR", ...
    "torqueRL", "rpmRL", "torqueRR", "rpmRR", "speed", "accelX"];
missing = strings(0,1);
for k = 1:numel(required)
    field = required(k);
    if isempty(c.(field)) || sum(isfinite(c.(field))) < 3
        missing(end+1) = field; %#ok<AGROW>
    end
end

% Fill short sensor gaps for event detection. Complete interpolation is done
% only after a session passes readiness checks and immediately before export.
for k = 1:numel(required)
    field = required(k);
    if ~isempty(c.(field))
        c.(field) = fillShortMissing(c.(field), t, cfg.MaxInterpolationGap_s);
    end
end
if ~isempty(c.drsState)
    c.drsState = round(fillShortMissing(c.drsState, t, cfg.MaxInterpolationGap_s));
end

session.channels = c;
session.estimatedChannels = unique(estimated, 'stable');
readiness = evaluateOptimizerReadiness(session, cfg);
end

function readiness = evaluateOptimizerReadiness(session, cfg)
c = session.channels;
required = ["tempFLRaw", "tempFRRaw", "tempRLRaw", "tempRRRaw", ...
    "frontPressureRaw", "rearPressureRaw", ...
    "torqueFL", "rpmFL", "torqueFR", "rpmFR", ...
    "torqueRL", "rpmRL", "torqueRR", "rpmRR", "speed", "accelX"];

missing = strings(0,1);
for k = 1:numel(required)
    field = required(k);
    if isempty(c.(field)) || sum(isfinite(c.(field))) < 3
        missing(end+1) = field; %#ok<AGROW>
    end
end

frontTemperatureSamples = max([finiteCount(c.tempFLRaw), finiteCount(c.tempFRRaw)]);
rearTemperatureSamples = max([finiteCount(c.tempRLRaw), finiteCount(c.tempRRRaw)]);

readiness = struct();
readiness.ready = isempty(missing);
readiness.missing = missing;
readiness.estimated = session.estimatedChannels;
readiness.frontTemperatureSamples = frontTemperatureSamples;
readiness.rearTemperatureSamples = rearTemperatureSamples;
readiness.hasTemperatureData = ...
    frontTemperatureSamples >= cfg.MinimumTemperatureSamples && ...
    rearTemperatureSamples >= cfg.MinimumTemperatureSamples;
end

function [session, info] = cleanTemperatureStartup(session, cfg)
% Clean EMI-related temperature spikes and trim to the first stable window.
info = struct( ...
    'stableStartFound', false, ...
    'trimmed', false, ...
    'trimmedSamples', 0, ...
    'trimmedDuration_s', 0, ...
    'originalStartTime_s', NaN, ...
    'newStartTime_s', NaN, ...
    'startTemperature_C', NaN, ...
    'validSensorCount', 0, ...
    'spikesReplaced', 0, ...
    'usedSpeedGate', false, ...
    'selectedSensors', strings(0,1), ...
    'rejectedSensors', strings(0,1), ...
    'reason', "disabled");

if ~cfg.TemperatureStartup.Enable
    return;
end

if isempty(session.time) || numel(session.time) < 3
    info.reason = "too_few_samples";
    return;
end

t = double(session.time(:));
dt = median(diff(t), 'omitnan');
if ~isfinite(dt) || dt <= 0
    info.reason = "invalid_sample_time";
    return;
end

temperatureFields = ["tempFLRaw", "tempFRRaw", "tempRLRaw", "tempRRRaw"];
tempRaw = nan(numel(t), numel(temperatureFields));
temperaturePresent = false(1,numel(temperatureFields));
for k = 1:numel(temperatureFields)
    value = session.channels.(temperatureFields(k));
    if ~isempty(value) && numel(value) == numel(t)
        tempRaw(:,k) = double(value(:));
        temperaturePresent(k) = true;
    end
end

tempC = 0.246 .* (tempRaw - 406);
plausible = tempC >= cfg.TemperatureStartup.MinimumPlausible_C & ...
    tempC <= cfg.TemperatureStartup.MaximumPlausible_C;
tempC(~plausible) = NaN;

[tempC, spikesReplaced] = filterTemperatureSpikes(tempC, t, cfg.TemperatureStartup);
info.spikesReplaced = spikesReplaced;

% Return cleaned temperature values to the raw ADC scale used by the optimizer.
cleanedRaw = tempC ./ 0.246 + 406;
for k = 1:numel(temperatureFields)
    if temperaturePresent(k)
        session.channels.(temperatureFields(k)) = cleanedRaw(:,k);
    end
end

windowSamples = max(3, round(cfg.TemperatureStartup.StabilityWindow_s / dt));
confirmationSamples = max(1, round(cfg.TemperatureStartup.ConfirmationWindow_s / dt));
searchEnd = find(t <= t(1) + cfg.TemperatureStartup.MaximumSearchTime_s, 1, 'last');
if isempty(searchEnd)
    searchEnd = numel(t);
end
searchEnd = min(searchEnd, numel(t));

if searchEnd < windowSamples
    info.reason = "search_window_too_short";
    return;
end

speed = [];
if ~isempty(session.channels.speed) && numel(session.channels.speed) == numel(t)
    speed = abs(double(session.channels.speed(:)));
end

% First try to locate a stable window while the car is stationary. If speed
% is unavailable, or no stationary stable window exists, retry without speed.
useSpeedOptions = false;
if ~isempty(speed) && sum(isfinite(speed)) >= windowSamples
    useSpeedOptions = [true, false];
end

stableIndex = [];
stableSensors = false(1,size(tempC,2));
usedSpeedGate = false;
for option = 1:numel(useSpeedOptions)
    useSpeedGate = useSpeedOptions(option);
    [stableIndex, stableSensors] = findStableTemperatureWindow( ...
        tempC, speed, windowSamples, confirmationSamples, searchEnd, dt, ...
        cfg.TemperatureStartup, useSpeedGate);
    if ~isempty(stableIndex)
        usedSpeedGate = useSpeedGate;
        break;
    end
end

if isempty(stableIndex)
    info.reason = "no_persistent_stable_window";
    return;
end

windowEnd = min(numel(t), stableIndex + windowSamples - 1);

% Reject non-selected corner channels that are clearly inconsistent with
% the coherent pair used to establish the startup temperature. This keeps a
% failed/noisy left sensor from blocking valid right-side temperature data.
[tempC, rejectedSensors] = rejectDisagreeingTemperatureChannels( ...
    tempC, stableIndex, windowEnd, stableSensors, ...
    cfg.TemperatureStartup.MaxInitialChannelDisagreement_C, ...
    cfg.TemperatureStartup.RejectDisagreeingChannels);

% Return the final cleaned/rejected temperature values to the raw ADC scale.
cleanedRaw = tempC ./ 0.246 + 406;
for k = 1:numel(temperatureFields)
    if temperaturePresent(k)
        session.channels.(temperatureFields(k)) = cleanedRaw(:,k);
    end
end

stableValues = tempC(stableIndex:windowEnd, stableSensors);
startTemperature = median(stableValues(:), 'omitnan');

sensorLabels = ["FL", "FR", "RL", "RR"];
info.stableStartFound = true;
info.trimmed = stableIndex > 1;
info.trimmedSamples = stableIndex - 1;
info.originalStartTime_s = t(1);
info.newStartTime_s = t(stableIndex);
info.trimmedDuration_s = t(stableIndex) - t(1);
info.startTemperature_C = startTemperature;
info.validSensorCount = sum(stableSensors);
info.usedSpeedGate = usedSpeedGate;
info.selectedSensors = sensorLabels(stableSensors).';
info.rejectedSensors = sensorLabels(rejectedSensors).';
info.reason = "stable_window_found";

if stableIndex > 1
    channelNames = fieldnames(session.channels);
    for k = 1:numel(channelNames)
        fieldName = channelNames{k};
        value = session.channels.(fieldName);
        if isnumeric(value) && size(value,1) == numel(t)
            session.channels.(fieldName) = value(stableIndex:end,:);
        end
    end
    session.time = t(stableIndex:end);
else
    session.time = t;
end

% Reset to a zero-based time vector after trimming. The source time range is
% retained in info for diagnostics, while all channels remain synchronized.
session.time = session.time - session.time(1);
session.channels.time = session.time;
end

function [stableIndex, stableSensors] = findStableTemperatureWindow( ...
        tempC, speed, windowSamples, confirmationSamples, searchEnd, dt, ...
        startupCfg, useSpeedGate)
% Find a stable coherent front/rear pair. The columns are:
%   1 = FL, 2 = FR, 3 = RL, 4 = RR
%
% A key detail is that stability is evaluated by candidate pairs rather than
% by all stable channels at once. A failed left-side sensor can therefore be
% ignored while a coherent right-side FR/RR pair remains usable.

stableIndex = [];
stableSensors = false(1,size(tempC,2));
lastStart = searchEnd - windowSamples - confirmationSamples + 2;
if lastStart < 1
    return;
end

candidateGroups = temperatureCandidateGroups(startupCfg, size(tempC,2));

for firstIndex = 1:lastStart
    windowIndex = firstIndex:(firstIndex + windowSamples - 1);
    confirmStart = windowIndex(end) + 1;
    confirmEnd = min(size(tempC,1), confirmStart + confirmationSamples - 1);
    confirmIndex = confirmStart:confirmEnd;

    windowTemps = tempC(windowIndex,:);
    validCounts = sum(isfinite(windowTemps),1);
    enoughData = validCounts >= max(3, ceil(0.8*windowSamples));

    ranges = max(windowTemps,[],1,'omitnan') - min(windowTemps,[],1,'omitnan');

    % Use the overall linear trend rather than point-to-point differences.
    trendRates = nan(1,size(windowTemps,2));
    relativeTime = (0:(windowSamples-1)).' .* dt;
    for sensor = 1:size(windowTemps,2)
        valid = isfinite(windowTemps(:,sensor));
        if sum(valid) >= 3
            fitCoefficients = polyfit(relativeTime(valid), ...
                windowTemps(valid,sensor), 1);
            trendRates(sensor) = abs(fitCoefficients(1));
        end
    end

    individuallyStable = enoughData & ...
        ranges <= startupCfg.MaxWindowRange_C & ...
        trendRates <= startupCfg.MaxMedianRate_C_per_s;

    selectedGroup = false(1,size(tempC,2));

    for groupNumber = 1:numel(candidateGroups)
        group = candidateGroups{groupNumber};
        group = group(group <= size(tempC,2));

        if numel(group) < startupCfg.MinimumValidSensors || ...
                any(~individuallyStable(group))
            continue;
        end

        groupMask = false(1,size(tempC,2));
        groupMask(group) = true;

        sensorMedians = median(windowTemps(:,group),1,'omitnan');
        if any(~isfinite(sensorMedians))
            continue;
        end

        % Only compare the sensors in this candidate group. Previously all
        % four stable sensors were compared, so large erroneous FL/RL values
        % caused a valid FR/RR pair to be rejected.
        if (max(sensorMedians) - min(sensorMedians)) > ...
                startupCfg.MaxSensorSpread_C
            continue;
        end

        if ~isempty(confirmIndex)
            confirmMedians = median(tempC(confirmIndex,group),1,'omitnan');
            confirmationShift = abs(confirmMedians - sensorMedians);
            if any(~isfinite(confirmMedians)) || ...
                    median(confirmationShift,'omitnan') > ...
                    startupCfg.MaxConfirmationShift_C
                continue;
            end
        end

        selectedGroup = groupMask;
        break;
    end

    if sum(selectedGroup) < startupCfg.MinimumValidSensors
        continue;
    end

    if useSpeedGate
        speedWindow = speed(windowIndex);
        if sum(isfinite(speedWindow)) < ceil(0.5*numel(speedWindow)) || ...
                median(speedWindow,'omitnan') > startupCfg.MaximumSpeed_mps
            continue;
        end
    end

    stableIndex = firstIndex;
    stableSensors = selectedGroup;
    return;
end
end

function groups = temperatureCandidateGroups(startupCfg, nSensors)
% Return front/rear pairs in the preferred order.
% [2 4] is FR/RR (right side), [1 3] is FL/RL (left side).

rightSide = [2 4];
leftSide = [1 3];

preferredSide = lower(string(startupCfg.PreferredSide));
if preferredSide == "left"
    groups = {leftSide};
    opposite = rightSide;
else
    groups = {rightSide};
    opposite = leftSide;
end

if startupCfg.AllowOppositeSideFallback
    groups{end+1} = opposite;
end

if startupCfg.AllowMixedSideFallback
    % Mixed-side fallbacks still require one front and one rear sensor.
    groups{end+1} = [1 4];
    groups{end+1} = [2 3];
end

% Remove groups that reference unavailable columns.
keep = true(size(groups));
for k = 1:numel(groups)
    keep(k) = all(groups{k} <= nSensors);
end
groups = groups(keep);
end

function [tempC, rejected] = rejectDisagreeingTemperatureChannels( ...
        tempC, stableIndex, windowEnd, selectedSensors, ...
        maximumDisagreement_C, enabled)
% Mark a non-selected corner as invalid when its startup median is far from
% the selected coherent pair. The complete channel is set to NaN so the bad
% values are not exported or accidentally used later.

rejected = false(1,size(tempC,2));
if ~enabled || ~any(selectedSensors)
    return;
end

selectedValues = tempC(stableIndex:windowEnd,selectedSensors);
referenceTemperature = median(selectedValues(:),'omitnan');
if ~isfinite(referenceTemperature)
    return;
end

for sensor = 1:size(tempC,2)
    if selectedSensors(sensor)
        continue;
    end

    sensorMedian = median(tempC(stableIndex:windowEnd,sensor),'omitnan');
    if ~isfinite(sensorMedian) || ...
            abs(sensorMedian-referenceTemperature) > maximumDisagreement_C
        tempC(:,sensor) = NaN;
        rejected(sensor) = true;
    end
end
end

function [filteredC, replacedCount] = filterTemperatureSpikes(tempC, time, startupCfg)
filteredC = tempC;
replacedCount = 0;
if ~startupCfg.SpikeFilter.Enable || isempty(tempC)
    return;
end

dt = median(diff(time),'omitnan');
if ~isfinite(dt) || dt <= 0
    return;
end

windowSamples = max(3, round(startupCfg.SpikeFilter.Window_s/dt));
if mod(windowSamples,2) == 0
    windowSamples = windowSamples + 1;
end

if startupCfg.SpikeFilter.ApplyToEntireFile
    filterEnd = size(tempC,1);
else
    filterEnd = find(time <= time(1) + startupCfg.MaximumSearchTime_s, 1, 'last');
    if isempty(filterEnd)
        filterEnd = size(tempC,1);
    end
end

for channel = 1:size(tempC,2)
    x = tempC(1:filterEnd,channel);
    if sum(isfinite(x)) < 3
        continue;
    end
    localMedian = movmedian(x, windowSamples, 'omitnan', 'Endpoints', 'shrink');
    deviation = abs(x-localMedian);
    localMad = movmedian(deviation, windowSamples, 'omitnan', 'Endpoints', 'shrink');
    threshold = max(startupCfg.SpikeFilter.MinimumDeviation_C, ...
        startupCfg.SpikeFilter.MadFactor .* 1.4826 .* localMad);
    spike = isfinite(x) & isfinite(localMedian) & deviation > threshold;
    x(spike) = localMedian(spike);
    filteredC(1:filterEnd,channel) = x;
    replacedCount = replacedCount + sum(spike);
end
end

function value = invalidateOutside(value, minimum, maximum)
if isempty(value)
    return;
end
value = double(value(:));
value(value < minimum | value > maximum) = NaN;
end

function [left, right, estimated] = completeAxlePair(left, right, leftName, rightName, estimated)
if isempty(left) && ~isempty(right)
    left = right;
    estimated(end+1) = leftName + "_from_" + rightName; %#ok<AGROW>
elseif isempty(right) && ~isempty(left)
    right = left;
    estimated(end+1) = rightName + "_from_" + leftName; %#ok<AGROW>
end
end

function matrix = collectAvailableChannels(channels, fields, nRows)
matrix = [];
for k = 1:numel(fields)
    value = channels.(fields(k));
    if ~isempty(value)
        value = double(value(:));
        if numel(value) == nRows
            matrix(:,end+1) = value; %#ok<AGROW>
        end
    end
end
end

function filled = fillShortMissing(value, time, maximumGap_s)
filled = double(value(:));
time = double(time(:));
missing = ~isfinite(filled);
if ~any(missing) || sum(~missing) < 2
    return;
end

runs = logicalRuns(missing);
for r = 1:size(runs,1)
    first = runs(r,1);
    last = runs(r,2);
    if first == 1 || last == numel(filled)
        continue;
    end
    gapDuration = time(last+1) - time(first-1);
    if gapDuration <= maximumGap_s
        idx = first:last;
        filled(idx) = interp1(time([first-1,last+1]), ...
            filled([first-1,last+1]), time(idx), 'linear');
    end
end
end

function filled = fillAllMissing(value, time)
filled = double(value(:));
time = double(time(:));
good = isfinite(filled) & isfinite(time);
if ~any(good)
    return;
elseif sum(good) == 1
    filled(:) = filled(find(good,1));
else
    filled = interp1(time(good), filled(good), time, 'linear', 'extrap');
end
end

%% Segment detection
function segments = detectSegments(session, cfg)
t = session.time(:);
c = session.channels;
if numel(t) < 3
    segments = struct([]);
    return;
end

dt = median(diff(t), 'omitnan');
if ~isfinite(dt) || dt <= 0
    segments = struct([]);
    return;
end

movementMask = buildMovementMask(c, t, cfg);

switch session.kind
    case "endurance_stint"
        settings = cfg.Endurance;
        mask = movementMask;
    case "autocross_run"
        settings = cfg.Autocross;
        mask = movementMask;
    case "acceleration_run"
        settings = cfg.Acceleration;
        mask = movementMask;
    case "braking_event"
        settings = cfg.Braking;
        mask = buildBrakingMask(c, t, cfg);
    otherwise
        settings = cfg.Moving;
        mask = movementMask;
end

% --- Noise robustness -----------------------------------------------
% Drop isolated blips (a handful of noisy samples) BEFORE bridging, so a
% cluster of brief sensor glitches can't be stitched together into
% something that merely spans long enough to pass MinimumDuration_s.
prebridgeMinSamples = max(2, round(cfg.NoiseRejection.MinRawBlip_s/dt));
mask = removeShortTrueRuns(mask, prebridgeMinSamples);
rawMaskBeforeBridge = mask;

mask = bridgeShortFalseRuns(mask, round(settings.MergeIdleGap_s/dt));

% Reject bridged runs that are mostly gap and only briefly touched real
% activity (e.g. a couple of seconds of a spiking sensor followed by a
% long idle gap, bridged together with another such spike). A genuine
% driving stint has activity sustained throughout its span, not just a
% low-duty-cycle sprinkling of it.
mask = enforceActivityDensity(mask, rawMaskBeforeBridge, cfg.NoiseRejection.MinActivityDensity);

mask = removeShortTrueRuns(mask, round(settings.MinimumDuration_s/dt));
runs = logicalRuns(mask);

if isempty(runs)
    segments = struct([]);
    return;
end

% Retain data before the detected driving period and after it so rotor
% cold-start and cooldown behavior remain available to BrakeCoeffOptimizer.
% Padding reaches back toward the END of the PREVIOUS stint (or the start
% of the file, for the first stint) and forward toward the START of the
% NEXT stint (or the end of the file, for the last stint), each bounded
% by a generous cap (settings.PrePadding_s / PostCooldown_s) so an
% unreasonably long idle gap between stints isn't pulled in wholesale.
% When two stints are close enough that their padded windows overlap,
% mergeOverlappingSegments below stitches them into one continuous file.
prePaddingCap = round(settings.PrePadding_s/dt);
postCooldownCap = round(settings.PostCooldown_s/dt);

nRuns = size(runs,1);
segments = repmat(struct('startIndex',0,'endIndex',0), nRuns, 1);
for r = 1:nRuns
    if r == 1
        prevRunEnd = 0;
    else
        prevRunEnd = runs(r-1,2);
    end
    if r == nRuns
        nextRunStart = numel(t)+1;
    else
        nextRunStart = runs(r+1,1);
    end

    segments(r).startIndex = max([1, prevRunEnd+1, runs(r,1)-prePaddingCap]);
    segments(r).endIndex   = min([numel(t), nextRunStart-1, runs(r,2)+postCooldownCap]);
end

% If a cooldown window reaches the following stint, retain the whole
% transition as one continuous segment rather than duplicating samples.
segments = mergeOverlappingSegments(segments);
end

function mask = enforceActivityDensity(bridgedMask, rawMask, minDensity)
% Drops any bridged run whose fraction of originally-active samples
% (before gap-bridging) falls below minDensity. This is what catches
% "a couple of seconds of spike, long gap, another couple of seconds of
% spike" patterns that would otherwise survive purely on bridged span
% length.
bridgedMask = logical(bridgedMask(:));
mask = bridgedMask;
if minDensity <= 0
    return;
end
rawMask = logical(rawMask(:));
runs = logicalRuns(bridgedMask);
for r = 1:size(runs,1)
    first = runs(r,1);
    last  = runs(r,2);
    density = sum(rawMask(first:last)) / (last-first+1);
    if density < minDensity
        mask(first:last) = false;
    end
end
end

function settings = getSegmentSettings(kind, cfg)
switch kind
    case "endurance_stint"
        settings = cfg.Endurance;
    case "autocross_run"
        settings = cfg.Autocross;
    case "acceleration_run"
        settings = cfg.Acceleration;
    case "braking_event"
        settings = cfg.Braking;
    otherwise
        settings = cfg.Moving;
end
end

function mask = buildMovementMask(c, t, cfg)
n = numel(t);
hasSpeed = ~isempty(c.speed) && sum(isfinite(c.speed) & c.speed > 0.1) >= max(10, round(0.005*n));

if hasSpeed
    mask = isfinite(c.speed) & c.speed >= cfg.MovingSpeedThreshold_mps;
else
    mask = false(n,1);
end

activity = false(n,1);
if ~isempty(c.accelX)
    accel = c.accelX(:);
    baseline = median(accel(isfinite(accel)), 'omitnan');
    accelFlag = isfinite(accel) & abs(accel-baseline) >= cfg.AccelActivityThreshold_mps2;
    windowSamples = max(1, round(cfg.ActivityWindow_s/median(diff(t), 'omitnan')));
    activity = movmean(double(accelFlag), windowSamples, 'Endpoints', 'shrink') ...
        >= cfg.ActivityFractionThreshold;
end

brakeActivity = false(n,1);
if ~isempty(c.frontPressureRaw)
    baseFront = sensorBaseline(c.frontPressureRaw);
    brakeActivity = brakeActivity | c.frontPressureRaw > baseFront + cfg.BrakeActivityDelta_ADC;
end
if ~isempty(c.rearPressureRaw)
    baseRear = sensorBaseline(c.rearPressureRaw);
    brakeActivity = brakeActivity | c.rearPressureRaw > baseRear + cfg.BrakeActivityDelta_ADC;
end
if any(brakeActivity)
    windowSamples = max(1, round(cfg.ActivityWindow_s/median(diff(t), 'omitnan')));
    brakeActivity = movmean(double(brakeActivity), windowSamples, 'Endpoints', 'shrink') > 0.1;
end

if hasSpeed
    % Activity helps bridge brief velocity dropouts without retaining long
    % stationary periods.
    mask = mask | (activity & brakeActivity);
else
    mask = activity | brakeActivity;
end
mask = logical(mask(:));
end

function mask = buildBrakingMask(c, t, cfg)
n = numel(t);
mask = false(n,1);

frontActive = false(n,1);
rearActive = false(n,1);
if ~isempty(c.frontPressureRaw)
    frontActive = c.frontPressureRaw > sensorBaseline(c.frontPressureRaw) + cfg.BrakeActivityDelta_ADC;
end
if ~isempty(c.rearPressureRaw)
    rearActive = c.rearPressureRaw > sensorBaseline(c.rearPressureRaw) + cfg.BrakeActivityDelta_ADC;
end
mask = frontActive | rearActive;

if ~isempty(c.accelX)
    mask = mask & isfinite(c.accelX) & c.accelX <= cfg.Braking.MinimumDecel_mps2;
end
if ~isempty(c.speed) && sum(isfinite(c.speed)) > 10
    mask = mask & isfinite(c.speed) & c.speed >= cfg.Braking.MinimumSpeed_mps;
end

windowSamples = max(1, round(0.2/median(diff(t), 'omitnan')));
mask = movmean(double(mask), windowSamples, 'Endpoints', 'shrink') > 0.25;
mask = logical(mask(:));
end

function baseline = sensorBaseline(value)
value = double(value(:));
valid = value(isfinite(value) & value > 100 & value < 2000);
if isempty(valid)
    valid = value(isfinite(value));
end
if isempty(valid)
    baseline = 0;
    return;
end
valid = sort(valid);
cutoff = max(1, round(0.40*numel(valid)));
baseline = median(valid(1:cutoff));
end

function mask = bridgeShortFalseRuns(mask, maximumGapSamples)
mask = logical(mask(:));
if maximumGapSamples <= 0
    return;
end
falseRuns = logicalRuns(~mask);
for r = 1:size(falseRuns,1)
    first = falseRuns(r,1);
    last = falseRuns(r,2);
    if first > 1 && last < numel(mask) && (last-first+1) <= maximumGapSamples
        mask(first:last) = true;
    end
end
end

function mask = removeShortTrueRuns(mask, minimumSamples)
mask = logical(mask(:));
minimumSamples = max(1, minimumSamples);
trueRuns = logicalRuns(mask);
for r = 1:size(trueRuns,1)
    first = trueRuns(r,1);
    last = trueRuns(r,2);
    if (last-first+1) < minimumSamples
        mask(first:last) = false;
    end
end
end

function runs = logicalRuns(mask)
mask = logical(mask(:));
changes = diff([false; mask; false]);
starts = find(changes == 1);
ends = find(changes == -1)-1;
runs = [starts, ends];
end

function segments = mergeOverlappingSegments(segments)
if numel(segments) <= 1
    return;
end
merged = segments(1);
for k = 2:numel(segments)
    if segments(k).startIndex <= merged(end).endIndex + 1
        merged(end).endIndex = max(merged(end).endIndex, segments(k).endIndex);
    else
        merged(end+1) = segments(k); %#ok<AGROW>
    end
end
segments = merged;
end

%% Optimizer formatting and writing
function [out, headers, formatName, info] = buildCombinedOptimizerMatrix(session, segments, cfg)
% Concatenate the detected segments while compressing the removed idle gaps.
% Each segment keeps its original internal sample timing; only the gap between
% segments is replaced with one nominal sample interval.

out = [];
headers = strings(1,0);
formatName = "";
sourceRanges = strings(numel(segments),1);
nominalDt = median(diff(session.time), 'omitnan');
if ~isfinite(nominalDt) || nominalDt <= 0
    nominalDt = 0.1;
end

for s = 1:numel(segments)
    seg = segments(s);
    [segmentMatrix, segmentHeaders, segmentFormat] = ...
        buildOptimizerMatrix(session, seg, cfg);

    if isempty(out)
        segmentMatrix(:,1) = segmentMatrix(:,1);
        headers = segmentHeaders;
        formatName = segmentFormat;
    else
        segmentMatrix(:,1) = segmentMatrix(:,1) + out(end,1) + nominalDt;
    end

    out = [out; segmentMatrix]; %#ok<AGROW>
    sourceRanges(s) = sprintf('%.3f-%.3f s', ...
        session.time(seg.startIndex), session.time(seg.endIndex));
end

info = struct();
info.firstSourceTime = session.time(segments(1).startIndex);
info.lastSourceTime = session.time(segments(end).endIndex);
info.outputDuration = out(end,1) - out(1,1);
info.sourceRanges = strjoin(sourceRanges, ', ');
end

function [out, headers, formatName] = buildOptimizerMatrix(session, segment, ~)
idx = segment.startIndex:segment.endIndex;
t = session.time(idx);
t = t-t(1);
c = session.channels;

fieldsToFill = ["tempFLRaw", "tempFRRaw", "tempRLRaw", "tempRRRaw", ...
    "frontPressureRaw", "rearPressureRaw", "steerPct", ...
    "accelPedalA", "accelPedalB", "torqueFL", "rpmFL", ...
    "torqueFR", "rpmFR", "torqueRL", "rpmRL", "torqueRR", "rpmRR", ...
    "speed", "accelX"];
for k = 1:numel(fieldsToFill)
    field = fieldsToFill(k);
    if isempty(c.(field))
        % Preserve the optimizer's fixed positional layout without inventing
        % a value for a channel that was not present in the source data.
        c.(field) = nan(size(session.time));
    else
        c.(field) = fillAllMissing(c.(field), session.time);
    end
end

hasDRS = ~isempty(c.drsState) && sum(isfinite(c.drsState)) >= 3;
reserved = zeros(numel(idx),1);

if hasDRS
    drs = round(fillAllMissing(c.drsState, session.time));
    out = [t, ...                                % 1
        c.tempFLRaw(idx), ...                    % 2
        c.tempFRRaw(idx), ...                    % 3
        c.tempRLRaw(idx), ...                    % 4
        c.tempRRRaw(idx), ...                    % 5
        drs(idx), ...                            % 6
        c.steerPct(idx), ...                     % 7
        c.frontPressureRaw(idx), ...             % 8
        c.accelPedalA(idx), ...                  % 9
        c.accelPedalB(idx), ...                  % 10
        reserved, ...                            % 11
        c.rearPressureRaw(idx), ...              % 12
        c.torqueFL(idx), c.rpmFL(idx), ...        % 13-14
        c.torqueFR(idx), c.rpmFR(idx), ...        % 15-16
        c.torqueRL(idx), c.rpmRL(idx), ...        % 17-18
        c.torqueRR(idx), c.rpmRR(idx), ...        % 19-20
        c.speed(idx), c.accelX(idx)];             % 21-22
    headers = ["Time_s", "FL_RotorTemp_RAW", "FR_RotorTemp_RAW", ...
        "RL_RotorTemp_RAW", "RR_RotorTemp_RAW", "DRS_State", ...
        "Steering_pct", "FrontBrakePressure_RAW", "AccelPedalA_pct", ...
        "AccelPedalB_pct", "Reserved", "RearBrakePressure_RAW", ...
        "FL_FeedbackTorque_pctMn", "FL_FeedbackVelocity_rpm", ...
        "FR_FeedbackTorque_pctMn", "FR_FeedbackVelocity_rpm", ...
        "RL_FeedbackTorque_pctMn", "RL_FeedbackVelocity_rpm", ...
        "RR_FeedbackTorque_pctMn", "RR_FeedbackVelocity_rpm", ...
        "Velocity_mps", "AccelX_mps2"];
    formatName = "drs_22_column";
else
    out = [t, ...                                % 1
        c.tempFLRaw(idx), ...                    % 2
        c.tempFRRaw(idx), ...                    % 3
        c.tempRLRaw(idx), ...                    % 4
        c.tempRRRaw(idx), ...                    % 5
        c.steerPct(idx), ...                     % 6
        c.frontPressureRaw(idx), ...             % 7
        c.accelPedalA(idx), ...                  % 8
        c.accelPedalB(idx), ...                  % 9
        reserved, ...                            % 10
        c.rearPressureRaw(idx), ...              % 11
        c.torqueFL(idx), c.rpmFL(idx), ...        % 12-13
        c.torqueFR(idx), c.rpmFR(idx), ...        % 14-15
        c.torqueRL(idx), c.rpmRL(idx), ...        % 16-17
        c.torqueRR(idx), c.rpmRR(idx), ...        % 18-19
        c.speed(idx), c.accelX(idx)];             % 20-21
    headers = ["Time_s", "FL_RotorTemp_RAW", "FR_RotorTemp_RAW", ...
        "RL_RotorTemp_RAW", "RR_RotorTemp_RAW", "Steering_pct", ...
        "FrontBrakePressure_RAW", "AccelPedalA_pct", "AccelPedalB_pct", ...
        "Reserved", "RearBrakePressure_RAW", ...
        "FL_FeedbackTorque_pctMn", "FL_FeedbackVelocity_rpm", ...
        "FR_FeedbackTorque_pctMn", "FR_FeedbackVelocity_rpm", ...
        "RL_FeedbackTorque_pctMn", "RL_FeedbackVelocity_rpm", ...
        "RR_FeedbackTorque_pctMn", "RR_FeedbackVelocity_rpm", ...
        "Velocity_mps", "AccelX_mps2"];
    formatName = "no_drs_21_column";
end

% NaN is intentionally retained for channels absent from the raw export.
% Convert any accidental infinities to NaN so the text file remains readable.
out(isinf(out)) = NaN;
end

function writeOptimizerTextFile(outputPath, matrix, headers, metadata)
fid = fopen(outputPath, 'w');
if fid < 0
    error('Could not open output file for writing: %s', outputPath);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# F34 BrakeDataPreprocessor optimizer-compatible export\n');
fprintf(fid, '# Source session: %s\n', metadata.sessionKey);
if isfield(metadata, 'segmentCount') && metadata.segmentCount > 1
    fprintf(fid, '# Segment type: %s | combined %d segments | source ranges: %s\n', ...
        metadata.kind, metadata.segmentCount, metadata.sourceRanges);
else
    fprintf(fid, '# Segment type: %s | source time %.3f to %.3f s\n', ...
        metadata.kind, metadata.startTime, metadata.endTime);
end
fprintf(fid, '# Source files: %s\n', strjoin(metadata.sourceFiles, ' | '));
if isempty(metadata.estimatedChannels)
    estimatedText = 'none';
else
    estimatedText = strjoin(metadata.estimatedChannels, ', ');
end
if ~isfield(metadata, 'missingChannels') || isempty(metadata.missingChannels)
    missingText = 'none';
else
    missingText = strjoin(metadata.missingChannels, ', ');
end
startupText = 'not-run';
if isfield(metadata, 'startupInfo')
    s = metadata.startupInfo;
    if s.stableStartFound
        startupText = sprintf('trimmed %.2f s, start %.1f C, spikes %d', ...
            s.trimmedDuration_s, s.startTemperature_C, s.spikesReplaced);
    else
        startupText = sprintf('no stable start (%s)', s.reason);
    end
end
fprintf(fid, ['# Format: %s | estimated: %s | missing-as-NaN: %s | ', ...
    'startup-cleaning: %s\n'], ...
    metadata.formatName, estimatedText, missingText, startupText);
fprintf(fid, '%s\n', strjoin(headers, sprintf('\t')));

nCols = size(matrix,2);
format = [repmat('%.10g\t', 1, nCols-1), '%.10g\n'];
for r = 1:size(matrix,1)
    fprintf(fid, format, matrix(r,:));
end
clear cleanup
end

function safe = sanitizeFileName(text)
safe = regexprep(string(text), '[^A-Za-z0-9_-]+', '_');
safe = regexprep(safe, '_+', '_');
safe = strip(safe, '_');
if strlength(safe) == 0
    safe = "session";
end
end

%% Manifest
function rows = appendManifestRow(rows, session, readiness, segmentNumber, ...
        startTime, endTime, duration, outputFile, status, segmentsCombined)
row = struct();
row.SessionKey = string(session.key);
row.SessionType = string(session.kind);
row.Status = string(status);
row.SegmentNumber = segmentNumber;
row.SegmentsCombined = segmentsCombined;
row.SourceStart_s = startTime;
row.SourceEnd_s = endTime;
row.Duration_s = duration;
row.OutputFile = string(outputFile);
row.MissingChannels = strjoin(readiness.missing, '; ');
row.EstimatedChannels = strjoin(readiness.estimated, '; ');
row.SourceFiles = strjoin(session.sourceFiles, '; ');

if isempty(rows)
    rows = row;
else
    rows(end+1) = row; %#ok<AGROW>
end
end

function manifest = makeManifestTable(rows)
if isempty(rows)
    manifest = table(strings(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
        zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), strings(0,1), ...
        strings(0,1), strings(0,1), strings(0,1), ...
        'VariableNames', {'SessionKey','SessionType','Status','SegmentNumber', ...
        'SegmentsCombined','SourceStart_s','SourceEnd_s','Duration_s','OutputFile', ...
        'MissingChannels','EstimatedChannels','SourceFiles'});
else
    manifest = struct2table(rows);
end
end
