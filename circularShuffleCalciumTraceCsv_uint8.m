%% Circularly shuffle calcium traces from a CSV and save a uint8 HDF5 cube.
% Required caller input:
%   traceCsvPath = path to a CSV with columns like cell_0, cell_1, ...
%
% Example:
%   traceCsvPath = '/path/to/K_Ca_traces_filtered_origHz.csv';
%   run('circularShuffleCalciumTraceCsv_uint8.m');
%
% Optional caller overrides:
%   numShuffles = number of circularly shifted shuffles to create (default: 1000)
%   cellPrefix = prefix used to identify trace columns (default: 'cell_')
%   outputDir = directory for HDF5/metadata outputs (default: input parent)
%   outputPrefix = output file stem prefix (default: '<inputBase>_circularShuffledTracesUint8')
%   maxChunkBytes = target bytes per write chunk (default: 128 * 1024 * 1024)
%   shuffleBatchSize = explicit shuffle batch size override
%   storageChunkShuffles = HDF5 chunk depth for shuffle axis (default: 1)
%   compressionLevel = HDF5 deflate level 0-9 (default: 9)
%   runTimestamp = output timestamp override (default: datestr(now, 'yyyymmdd_HHMMSS'))

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
    error('Non-finite values were found in the selected cell columns. Clean the input CSV before shuffling.');
end

[numFrames, numCells] = size(traceMatrix);
if numFrames < 2
    error('At least 2 rows are required to perform circular shuffles.');
end

[inputDir, inputBase, ~] = fileparts(traceCsvPath);
if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = inputDir;
end

if ~exist('outputPrefix', 'var') || isempty(outputPrefix)
    outputPrefix = [inputBase '_circularShuffledTracesUint8'];
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

outputStem = [outputPrefix '_' runTimestamp];
h5FilePath = fullfile(outputDir, [outputStem '.h5']);
cellMetadataPath = fullfile(outputDir, [outputStem '_cellMetadata.csv']);
configPath = fullfile(outputDir, [outputStem '_config.json']);

traceMinPerCell = min(traceMatrix, [], 1);
traceMaxPerCell = max(traceMatrix, [], 1);
traceRangePerCell = traceMaxPerCell - traceMinPerCell;
zeroRangeCellMask = traceRangePerCell == 0;
traceRangeForQuantization = traceRangePerCell;
traceRangeForQuantization(zeroRangeCellMask) = 1;

quantizedTraceMatrix = quantizeTraceMatrixUint8(traceMatrix, traceMinPerCell, traceRangeForQuantization);

bytesPerShuffle = numFrames * numCells;
if ~exist('shuffleBatchSize', 'var') || isempty(shuffleBatchSize)
    shuffleBatchSize = max(1, floor(double(maxChunkBytes) / double(bytesPerShuffle)));
end
shuffleBatchSize = min(numShuffles, max(1, round(shuffleBatchSize)));
storageChunkShuffles = min([numShuffles, shuffleBatchSize, max(1, round(storageChunkShuffles))]);

recreateFileIfNeeded(h5FilePath);

h5create(h5FilePath, '/quantizedTraceMatrixUint8', [numFrames, numCells], 'Datatype', 'uint8');
h5create(h5FilePath, '/shuffledTracesUint8', [numFrames, numCells, Inf], ...
    'Datatype', 'uint8', ...
    'ChunkSize', [numFrames, numCells, storageChunkShuffles], ...
    'Deflate', compressionLevel);
h5create(h5FilePath, '/circularShiftByCellAndShuffle', [numCells, numShuffles], ...
    'Datatype', 'uint32', ...
    'ChunkSize', [numCells, min(numShuffles, shuffleBatchSize)], ...
    'Deflate', compressionLevel);
h5create(h5FilePath, '/traceMinPerCell', size(traceMinPerCell), 'Datatype', 'double');
h5create(h5FilePath, '/traceRangePerCell', size(traceRangePerCell), 'Datatype', 'double');
h5create(h5FilePath, '/zeroRangeCellMask', size(zeroRangeCellMask), 'Datatype', 'uint8');
h5create(h5FilePath, '/matrixShape', [2, 1], 'Datatype', 'uint64');

h5write(h5FilePath, '/quantizedTraceMatrixUint8', quantizedTraceMatrix);
h5write(h5FilePath, '/traceMinPerCell', traceMinPerCell);
h5write(h5FilePath, '/traceRangePerCell', traceRangePerCell);
h5write(h5FilePath, '/zeroRangeCellMask', uint8(zeroRangeCellMask));
h5write(h5FilePath, '/matrixShape', uint64([numFrames; numCells]));

for startShuffle = 1:shuffleBatchSize:numShuffles
    endShuffle = min(startShuffle + shuffleBatchSize - 1, numShuffles);
    currentChunkSize = endShuffle - startShuffle + 1;

    disp(sprintf('writing shuffles %d-%d of %d', startShuffle, endShuffle, numShuffles)); %#ok<DSPS>

    shuffledChunk = zeros(numFrames, numCells, currentChunkSize, 'uint8');
    shiftChunk = zeros(numCells, currentChunkSize, 'uint32');

    for shuffleIdx = 1:currentChunkSize
        shiftsThisShuffle = randi([0, numFrames - 1], numCells, 1);
        shiftChunk(:, shuffleIdx) = uint32(shiftsThisShuffle);

        for cellIdx = 1:numCells
            shuffledChunk(:, cellIdx, shuffleIdx) = circshift(quantizedTraceMatrix(:, cellIdx), shiftsThisShuffle(cellIdx));
        end
    end

    h5write(h5FilePath, '/shuffledTracesUint8', shuffledChunk, [1, 1, startShuffle], [numFrames, numCells, currentChunkSize]);
    h5write(h5FilePath, '/circularShiftByCellAndShuffle', shiftChunk, [1, startShuffle], [numCells, currentChunkSize]);
end

cellMetadataTable = table( ...
    (1:numCells)', ...
    string(cellColumns(:)), ...
    traceMinPerCell(:), ...
    traceMaxPerCell(:), ...
    traceRangePerCell(:), ...
    zeroRangeCellMask(:), ...
    'VariableNames', {'cellIndex', 'cellName', 'traceMin', 'traceMax', 'traceRange', 'zeroRange'});
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
configStruct.quantizationRule = 'per-cell min/max mapped to uint8 with dequantization value = (uint8 / 255) * traceRange + traceMin';
writeJsonConfig(configPath, configStruct);

disp('saved circular shuffle cube');
disp(h5FilePath);

function quantizedMatrix = quantizeTraceMatrixUint8(traceMatrix, traceMinPerCell, traceRangePerCell)
    [numFrames, numCells] = size(traceMatrix);
    quantizedMatrix = zeros(numFrames, numCells, 'uint8');
    for cellIdx = 1:numCells
        normalizedTrace = (traceMatrix(:, cellIdx) - traceMinPerCell(cellIdx)) ./ traceRangePerCell(cellIdx);
        normalizedTrace = max(0, min(1, normalizedTrace));
        quantizedMatrix(:, cellIdx) = uint8(round(255 * normalizedTrace));
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
