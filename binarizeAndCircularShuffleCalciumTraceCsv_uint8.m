%% Detect conservative calcium peaks, binarize them, and save shuffled uint8 cubes.
% Required caller input:
%   traceCsvPath = path to a CSV with columns like cell_0, cell_1, ...
%
% Example:
%   traceCsvPath = '/path/to/K_Ca_traces_filtered_origHz.csv';
%   run('binarizeAndCircularShuffleCalciumTraceCsv_uint8.m');
%
% Optional caller overrides:
%   numShuffles = number of circularly shifted shuffles to create (default: 1000)
%   cellPrefix = prefix used to identify trace columns (default: 'cell_')
%   outputDir = directory for HDF5/metadata outputs (default: input parent)
%   outputPrefix = output file stem prefix (default: '<inputBase>_binarizedPeakCircularShufflesUint8')
%   maxChunkBytes = target bytes per write chunk (default: 128 * 1024 * 1024)
%   shuffleBatchSize = explicit shuffle batch size override
%   storageChunkShuffles = HDF5 chunk depth for shuffle axis (default: 1)
%   compressionLevel = HDF5 deflate level 0-9 (default: 9)
%   runTimestamp = output timestamp override (default: datestr(now, 'yyyymmdd_HHMMSS'))
%   binaryMarkMode = 'peak' (default), 'onset', or 'event_window'
%   threshold = 2.5
%   smoothWindowSamples = 3
%   minPeakProminence = 1.5
%   minPeakDistanceSamples = 20
%   minPeakWidthSamples = 2
%   returnToBaselineTol = 0.4
%   returnToBaselineWindowSamples = 5
%   tol = 0.25
%   windowLenOnset = 3
%   windowLenOffset = 3
%   progressEveryCells = status update cadence while binarizing (default: 25)

if ~exist('traceCsvPath', 'var') || isempty(traceCsvPath)
    error('Set traceCsvPath to the calcium trace CSV before running this script.');
end

if ~exist('numShuffles', 'var') || isempty(numShuffles)
    numShuffles = 1000;
end

if ~exist('cellPrefix', 'var') || isempty(cellPrefix)
    cellPrefix = 'cell_';
end

if ~exist('maxChunkBytes', 'var') || isempty(maxChunkBytes)
    maxChunkBytes = 128 * 1024 * 1024;
end

if ~exist('storageChunkShuffles', 'var') || isempty(storageChunkShuffles)
    storageChunkShuffles = 1;
end

if ~exist('compressionLevel', 'var') || isempty(compressionLevel)
    compressionLevel = 9;
end

if ~exist('runTimestamp', 'var') || isempty(runTimestamp)
    runTimestamp = datestr(now, 'yyyymmdd_HHMMSS');
end

if ~exist('binaryMarkMode', 'var') || isempty(binaryMarkMode)
    binaryMarkMode = 'peak';
end

if ~exist('threshold', 'var') || isempty(threshold)
    threshold = 2.5;
end

if ~exist('smoothWindowSamples', 'var') || isempty(smoothWindowSamples)
    smoothWindowSamples = 3;
end

if ~exist('minPeakProminence', 'var') || isempty(minPeakProminence)
    minPeakProminence = 1.5;
end

if ~exist('minPeakDistanceSamples', 'var') || isempty(minPeakDistanceSamples)
    minPeakDistanceSamples = 20;
end

if ~exist('minPeakWidthSamples', 'var') || isempty(minPeakWidthSamples)
    minPeakWidthSamples = 2;
end

if ~exist('returnToBaselineTol', 'var') || isempty(returnToBaselineTol)
    returnToBaselineTol = 0.4;
end

if ~exist('returnToBaselineWindowSamples', 'var') || isempty(returnToBaselineWindowSamples)
    returnToBaselineWindowSamples = 5;
end

if ~exist('tol', 'var') || isempty(tol)
    tol = 0.25;
end

if ~exist('windowLenOnset', 'var') || isempty(windowLenOnset)
    windowLenOnset = 3;
end

if ~exist('windowLenOffset', 'var') || isempty(windowLenOffset)
    windowLenOffset = 3;
end

if ~exist('progressEveryCells', 'var') || isempty(progressEveryCells)
    progressEveryCells = 25;
end

validMarkModes = {'peak', 'onset', 'event_window'};
if ~ismember(binaryMarkMode, validMarkModes)
    error('binaryMarkMode must be one of: %s', strjoin(validMarkModes, ', '));
end

disp('loading trace CSV');
disp(traceCsvPath);
traceTable = readtable(traceCsvPath, 'VariableNamingRule', 'preserve');

variableNames = traceTable.Properties.VariableNames;
cellColumns = variableNames(startsWith(variableNames, cellPrefix));

if isempty(cellColumns)
    error('No columns beginning with "%s" were found in %s.', cellPrefix, traceCsvPath);
end

traceMatrix = table2array(traceTable(:, cellColumns));
traceMatrix = double(traceMatrix);

if any(~isfinite(traceMatrix), 'all')
    error('Non-finite values were found in the selected cell columns. Clean the input CSV before binarizing.');
end

[numFrames, numCells] = size(traceMatrix);
if numFrames < 2
    error('At least 2 rows are required to perform binarization and shuffling.');
end

[inputDir, inputBase, ~] = fileparts(traceCsvPath);
if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = inputDir;
end

if ~exist('outputPrefix', 'var') || isempty(outputPrefix)
    outputPrefix = [inputBase '_binarizedPeakCircularShufflesUint8'];
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

outputStem = [outputPrefix '_' runTimestamp];
h5FilePath = fullfile(outputDir, [outputStem '.h5']);
cellMetadataPath = fullfile(outputDir, [outputStem '_cellMetadata.csv']);
configPath = fullfile(outputDir, [outputStem '_config.json']);

detectOpts = struct();
detectOpts.threshold = threshold;
detectOpts.smoothWindowSamples = smoothWindowSamples;
detectOpts.minPeakProminence = minPeakProminence;
detectOpts.minPeakDistanceSamples = minPeakDistanceSamples;
detectOpts.minPeakWidthSamples = minPeakWidthSamples;
detectOpts.returnToBaselineTol = returnToBaselineTol;
detectOpts.returnToBaselineWindowSamples = returnToBaselineWindowSamples;
detectOpts.tol = tol;
detectOpts.windowLenOnset = windowLenOnset;
detectOpts.windowLenOffset = windowLenOffset;
detectOpts.binaryMarkMode = binaryMarkMode;

binarizedPeakMatrix = zeros(numFrames, numCells, 'uint8');
peakCountPerCell = zeros(numCells, 1);

for cellIdx = 1:numCells
    if mod(cellIdx, progressEveryCells) == 0 || cellIdx == 1 || cellIdx == numCells
        disp(sprintf('binarizing cell %d of %d', cellIdx, numCells)); %#ok<DSPS>
    end

    traceVector = traceMatrix(:, cellIdx);
    peakIdxs = detectPeakIndicesFromTrace(traceVector, detectOpts);
    peakCountPerCell(cellIdx) = numel(peakIdxs);

    if isempty(peakIdxs)
        continue;
    end

    switch binaryMarkMode
        case 'peak'
            binarizedPeakMatrix(peakIdxs, cellIdx) = 1;
        case 'onset'
            for peakNumber = 1:numel(peakIdxs)
                [onsetIdx, ~] = getPeakBoundsForTrace(peakIdxs(peakNumber), traceVector, traceVector, detectOpts);
                binarizedPeakMatrix(onsetIdx, cellIdx) = 1;
            end
        case 'event_window'
            for peakNumber = 1:numel(peakIdxs)
                [onsetIdx, offsetIdx] = getPeakBoundsForTrace(peakIdxs(peakNumber), traceVector, traceVector, detectOpts);
                binarizedPeakMatrix(onsetIdx:offsetIdx, cellIdx) = 1;
            end
    end
end

bytesPerShuffle = numFrames * numCells;
if ~exist('shuffleBatchSize', 'var') || isempty(shuffleBatchSize)
    shuffleBatchSize = max(1, floor(double(maxChunkBytes) / double(bytesPerShuffle)));
end
shuffleBatchSize = min(numShuffles, max(1, round(shuffleBatchSize)));
storageChunkShuffles = min([numShuffles, shuffleBatchSize, max(1, round(storageChunkShuffles))]);

recreateFileIfNeeded(h5FilePath);

h5create(h5FilePath, '/binarizedPeaks', [numFrames, numCells], 'Datatype', 'uint8');
h5create(h5FilePath, '/binarizedPeaksShuffled', [numFrames, numCells, Inf], ...
    'Datatype', 'uint8', ...
    'ChunkSize', [numFrames, numCells, storageChunkShuffles], ...
    'Deflate', compressionLevel);
h5create(h5FilePath, '/circularShiftByCellAndShuffle', [numCells, numShuffles], ...
    'Datatype', 'uint32', ...
    'ChunkSize', [numCells, min(numShuffles, shuffleBatchSize)], ...
    'Deflate', compressionLevel);
h5create(h5FilePath, '/peakCountPerCell', size(peakCountPerCell), 'Datatype', 'double');
h5create(h5FilePath, '/matrixShape', [2, 1], 'Datatype', 'uint64');

h5write(h5FilePath, '/binarizedPeaks', binarizedPeakMatrix);
h5write(h5FilePath, '/peakCountPerCell', peakCountPerCell);
h5write(h5FilePath, '/matrixShape', uint64([numFrames; numCells]));

for startShuffle = 1:shuffleBatchSize:numShuffles
    endShuffle = min(startShuffle + shuffleBatchSize - 1, numShuffles);
    currentChunkSize = endShuffle - startShuffle + 1;

    disp(sprintf('writing shuffled binary cube %d-%d of %d', startShuffle, endShuffle, numShuffles)); %#ok<DSPS>

    shuffledChunk = zeros(numFrames, numCells, currentChunkSize, 'uint8');
    shiftChunk = zeros(numCells, currentChunkSize, 'uint32');

    for shuffleIdx = 1:currentChunkSize
        shiftsThisShuffle = randi([0, numFrames - 1], numCells, 1);
        shiftChunk(:, shuffleIdx) = uint32(shiftsThisShuffle);

        for cellIdx = 1:numCells
            shuffledChunk(:, cellIdx, shuffleIdx) = circshift(binarizedPeakMatrix(:, cellIdx), shiftsThisShuffle(cellIdx));
        end
    end

    h5write(h5FilePath, '/binarizedPeaksShuffled', shuffledChunk, [1, 1, startShuffle], [numFrames, numCells, currentChunkSize]);
    h5write(h5FilePath, '/circularShiftByCellAndShuffle', shiftChunk, [1, startShuffle], [numCells, currentChunkSize]);
end

cellMetadataTable = table( ...
    (1:numCells)', ...
    string(cellColumns(:)), ...
    peakCountPerCell(:), ...
    'VariableNames', {'cellIndex', 'cellName', 'peakCount'});
writetable(cellMetadataTable, cellMetadataPath);

configStruct = struct();
configStruct.traceCsvPath = traceCsvPath;
configStruct.outputH5 = h5FilePath;
configStruct.outputCellMetadataCsv = cellMetadataPath;
configStruct.numFrames = numFrames;
configStruct.numCells = numCells;
configStruct.numShuffles = numShuffles;
configStruct.cellPrefix = cellPrefix;
configStruct.outputPrefix = outputPrefix;
configStruct.runTimestamp = runTimestamp;
configStruct.maxChunkBytes = maxChunkBytes;
configStruct.shuffleBatchSize = shuffleBatchSize;
configStruct.storageChunkShuffles = storageChunkShuffles;
configStruct.compressionLevel = compressionLevel;
configStruct.binaryMarkMode = binaryMarkMode;
configStruct.threshold = threshold;
configStruct.smoothWindowSamples = smoothWindowSamples;
configStruct.minPeakProminence = minPeakProminence;
configStruct.minPeakDistanceSamples = minPeakDistanceSamples;
configStruct.minPeakWidthSamples = minPeakWidthSamples;
configStruct.returnToBaselineTol = returnToBaselineTol;
configStruct.returnToBaselineWindowSamples = returnToBaselineWindowSamples;
configStruct.tol = tol;
configStruct.windowLenOnset = windowLenOnset;
configStruct.windowLenOffset = windowLenOffset;
writeJsonConfig(configPath, configStruct);

disp('saved binarized peak shuffle cube');
disp(h5FilePath);

function peakIdxs = detectPeakIndicesFromTrace(traceVector, detectOpts)
    traceVector = double(traceVector(:));
    workTrace = traceVector;

    if ~isempty(detectOpts.smoothWindowSamples) && detectOpts.smoothWindowSamples > 1
        workTrace = movmean(workTrace, detectOpts.smoothWindowSamples, 'Endpoints', 'shrink');
    end

    peakArgs = { ...
        'MinPeakHeight', detectOpts.threshold, ...
        'MinPeakDistance', max(1, round(detectOpts.minPeakDistanceSamples))};

    if ~isempty(detectOpts.minPeakProminence)
        peakArgs = [peakArgs, {'MinPeakProminence', detectOpts.minPeakProminence}]; %#ok<AGROW>
    end
    if ~isempty(detectOpts.minPeakWidthSamples)
        peakArgs = [peakArgs, {'MinPeakWidth', detectOpts.minPeakWidthSamples}]; %#ok<AGROW>
    end

    [~, peakIdxs] = findpeaks(workTrace, peakArgs{:});
    peakIdxs = peakIdxs(:);

    if isempty(peakIdxs) || isempty(detectOpts.returnToBaselineTol) || numel(peakIdxs) <= 1
        return;
    end

    acceptedPeaks = peakIdxs(1);
    for peakNumber = 2:numel(peakIdxs)
        candidateIdx = peakIdxs(peakNumber);
        if hasBaselineReturn(workTrace, acceptedPeaks(end), candidateIdx, detectOpts.returnToBaselineTol, detectOpts.returnToBaselineWindowSamples)
            acceptedPeaks(end + 1, 1) = candidateIdx; %#ok<AGROW>
        end
    end

    peakIdxs = acceptedPeaks;
end

function tf = hasBaselineReturn(workTrace, startIdx, stopIdx, tolerance, windowLength)
    if stopIdx <= startIdx + 1
        tf = true;
        return;
    end

    segment = abs(workTrace(startIdx + 1:stopIdx - 1));
    if isempty(segment)
        tf = true;
        return;
    end

    windowLength = max(1, round(windowLength));
    if numel(segment) < windowLength
        tf = all(segment <= tolerance);
        return;
    end

    baselineMask = segment <= tolerance;
    baselineRunLengths = conv(double(baselineMask), ones(windowLength, 1), 'valid');
    tf = any(baselineRunLengths == windowLength);
end

function [onsetIdx, offsetIdx] = getPeakBoundsForTrace(peakIdx, zTrace, rawTrace, detectOpts)
    zTrace = double(zTrace(:));
    rawTrace = double(rawTrace(:));
    numSamples = numel(zTrace);

    onsetIdx = 1;
    for idx = peakIdx:-1:detectOpts.windowLenOnset
        startIdx = idx - detectOpts.windowLenOnset + 1;
        if all(abs(zTrace(startIdx:idx)) < detectOpts.tol)
            onsetIdx = startIdx;
            break;
        end
    end

    baselineValue = rawTrace(onsetIdx);
    offsetIdx = numSamples;
    lastStartIdx = numSamples - detectOpts.windowLenOffset + 1;
    for idx = peakIdx:lastStartIdx
        windowIdx = idx:(idx + detectOpts.windowLenOffset - 1);
        if all(abs(zTrace(windowIdx)) < detectOpts.tol) && all(rawTrace(windowIdx) < baselineValue)
            offsetIdx = windowIdx(end);
            break;
        end
    end
end

function recreateFileIfNeeded(filePath)
    if isfile(filePath)
        delete(filePath);
    end
end

function writeJsonConfig(configPath, configStruct)
    fileId = fopen(configPath, 'w');
    if fileId == -1
        error('Could not open %s for writing.', configPath);
    end
    cleaner = onCleanup(@() fclose(fileId));
    fprintf(fileId, '%s\n', jsonencode(configStruct, PrettyPrint=true));
    clear cleaner;
end
