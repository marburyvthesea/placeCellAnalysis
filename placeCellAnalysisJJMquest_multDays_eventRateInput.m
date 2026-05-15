%% Split a combined multi-session event-rate dataframe and run place cell analysis per session.
% This variant starts from precomputed event-rate matrices such as
% instantaneousEventRate.csv or instantaneousEventRateOnsets.csv.
% It skips the signal-peaks and sliding-window event-rate steps and goes
% straight into the mutual-information workflow.
%
% Expected input:
%   alignedFile = path to a combined event-rate CSV
%
% Optional caller overrides:
%   metadataFile = companion GCAMP_with_velocity.csv path
%                  (default: sibling file next to alignedFile)
%   numShuffles = number of shuffle iterations (default: 1000)
%   numBins = number of spatial bins (default: 32)
%   binsubset = subset of bins used for MI subset output (default: 2:31)
%   maxChunkElements = target maximum elements per shuffle-results chunk (default: 5e8)
%   maxShufflesPerWrite = maximum shuffle count to accumulate before writing MI results (default: 100)
%   writeSplitSessionTables = whether to write one CSV per session (default: true)
%   sessionIdColumn = column used to split sessions (default: 'ezTrackOutput')
%   frameColumn = column used to sort rows within session (default: 'closestBehavCamFrameIdx')
%   xPositionColumn = column used for spatial position (default: 'X_coor')
%   runTimestamp = output timestamp override (default: datestr(now, 'yyyymmdd_HHMMSS'))

if ~exist('alignedFile', 'var') || isempty(alignedFile)
    error('Set alignedFile to the combined event-rate CSV path before running this script.');
end

if ~exist('numShuffles', 'var') || isempty(numShuffles)
    numShuffles = 1000;
end

if ~exist('numBins', 'var') || isempty(numBins)
    numBins = 32;
end

if ~exist('binsubset', 'var') || isempty(binsubset)
    binsubset = 2:31;
end

if ~exist('maxChunkElements', 'var') || isempty(maxChunkElements)
    maxChunkElements = 500000000;
end

if ~exist('maxShufflesPerWrite', 'var') || isempty(maxShufflesPerWrite)
    maxShufflesPerWrite = 100;
end

if ~exist('writeSplitSessionTables', 'var') || isempty(writeSplitSessionTables)
    writeSplitSessionTables = true;
end

if ~exist('sessionIdColumn', 'var') || isempty(sessionIdColumn)
    sessionIdColumn = 'ezTrackOutput';
end

if ~exist('frameColumn', 'var') || isempty(frameColumn)
    frameColumn = 'closestBehavCamFrameIdx';
end

if ~exist('xPositionColumn', 'var') || isempty(xPositionColumn)
    xPositionColumn = 'X_coor';
end

if ~exist('runTimestamp', 'var') || isempty(runTimestamp)
    runTimestamp = datestr(now, 'yyyymmdd_HHMMSS');
end

[inputDir, inputBase, ~] = fileparts(alignedFile);
if ~exist('metadataFile', 'var') || isempty(metadataFile)
    metadataFile = fullfile(inputDir, 'GCAMP_with_velocity.csv');
end

disp('loading event-rate input');
disp(alignedFile);
eventRateTable = readtable(alignedFile, 'VariableNamingRule', 'preserve');

eventRateVariableNames = eventRateTable.Properties.VariableNames;
cellColumns = eventRateVariableNames(startsWith(eventRateVariableNames, 'cell_'));

if isempty(cellColumns)
    error('No cell_* columns were found in %s.', alignedFile);
end

requiredColumns = {sessionIdColumn, xPositionColumn};
hasEmbeddedMetadata = all(ismember(requiredColumns, eventRateVariableNames));

if hasEmbeddedMetadata
    dataTable = eventRateTable;
else
    if ~isfile(metadataFile)
        error(['Event-rate input %s does not include %s/%s, and companion metadata file ' ...
            '%s was not found.'], alignedFile, sessionIdColumn, xPositionColumn, metadataFile);
    end

    disp('loading companion metadata');
    disp(metadataFile);
    metadataTable = readtable(metadataFile, 'VariableNamingRule', 'preserve');

    if height(eventRateTable) ~= height(metadataTable)
        error(['Row count mismatch between event-rate input (%d rows) and metadata file (%d rows). ' ...
            'These files must align frame-by-frame.'], height(eventRateTable), height(metadataTable));
    end

    metadataVariableNames = metadataTable.Properties.VariableNames;
    for colIdx = 1:numel(requiredColumns)
        if ~ismember(requiredColumns{colIdx}, metadataVariableNames)
            error('Required metadata column "%s" was not found in %s.', requiredColumns{colIdx}, metadataFile);
        end
    end

    metadataKeepMask = ~startsWith(metadataVariableNames, 'cell_');
    dataTable = [metadataTable(:, metadataKeepMask), eventRateTable(:, cellColumns)];
end

variableNames = dataTable.Properties.VariableNames;
for colIdx = 1:numel(requiredColumns)
    if ~ismember(requiredColumns{colIdx}, variableNames)
        error('Required column "%s" was not found after preparing %s.', requiredColumns{colIdx}, alignedFile);
    end
end

sessionIds = string(dataTable.(sessionIdColumn));
validRows = ~ismissing(sessionIds) & strlength(strtrim(sessionIds)) > 0;

if ~all(validRows)
    warning('Dropping %d rows with missing session identifiers.', sum(~validRows));
    dataTable = dataTable(validRows, :);
    sessionIds = sessionIds(validRows);
end

outputDir = fullfile(inputDir, [inputBase '_placeCellAnalysis_multDays_eventRateInput_' runTimestamp]);
splitDir = fullfile(outputDir, 'splitSessions');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

if writeSplitSessionTables && ~isfolder(splitDir)
    mkdir(splitDir);
end

sessionSources = unique(sessionIds, 'stable');
numSessions = numel(sessionSources);

manifest = table( ...
    zeros(numSessions, 1), ...
    strings(numSessions, 1), ...
    zeros(numSessions, 1), ...
    strings(numSessions, 1), ...
    strings(numSessions, 1), ...
    strings(numSessions, 1), ...
    strings(numSessions, 1), ...
    'VariableNames', {'sessionIndex', 'sessionSource', 'numRows', 'sessionStem', 'splitCsv', 'status', 'message'});

opts = struct();
opts.numShuffles = numShuffles;
opts.numBins = numBins;
opts.binsubset = binsubset;
opts.maxChunkElements = maxChunkElements;
opts.maxShufflesPerWrite = maxShufflesPerWrite;
opts.xPositionColumn = xPositionColumn;

for sessionIdx = 1:numSessions
    sessionSource = sessionSources(sessionIdx);
    sessionMask = sessionIds == sessionSource;
    sessionTable = dataTable(sessionMask, :);

    if ismember(frameColumn, sessionTable.Properties.VariableNames)
        sessionTable = sortrows(sessionTable, frameColumn);
    end

    sessionStem = buildSessionStem(inputBase, sessionIdx, sessionSource, runTimestamp);
    splitCsvPath = "";

    manifest.sessionIndex(sessionIdx) = sessionIdx;
    manifest.sessionSource(sessionIdx) = sessionSource;
    manifest.numRows(sessionIdx) = height(sessionTable);
    manifest.sessionStem(sessionIdx) = sessionStem;

    disp(['processing session ' num2str(sessionIdx) ' of ' num2str(numSessions)]);
    disp(sessionSource);

    try
        if writeSplitSessionTables
            splitCsvPath = fullfile(splitDir, [sessionStem '.csv']);
            writetable(sessionTable, splitCsvPath);
        end

        analyzeSingleSession(sessionTable, cellColumns, outputDir, sessionStem, opts);

        manifest.splitCsv(sessionIdx) = splitCsvPath;
        manifest.status(sessionIdx) = "ok";
        manifest.message(sessionIdx) = "";
    catch ME
        manifest.splitCsv(sessionIdx) = splitCsvPath;
        manifest.status(sessionIdx) = "error";
        manifest.message(sessionIdx) = string(ME.message);
        warning('Session %d failed: %s', sessionIdx, getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

manifestPath = fullfile(outputDir, [inputBase '__run' runTimestamp '_session_manifest.csv']);
writetable(manifest, manifestPath);

disp('multi-session event-rate place cell analysis completed');
disp(outputDir);

function analyzeSingleSession(sessionTable, cellColumns, outputDir, outputStem, opts)
    eventRateArray = table2array(sessionTable(:, cellColumns));
    x_position = sessionTable.(opts.xPositionColumn);

    [numFrames, numNeurons] = size(eventRateArray);
    if numFrames < 2
        error('Session %s has fewer than 2 frames.', outputStem);
    end

    if numNeurons < 1
        error('Session %s has no cell columns.', outputStem);
    end

    invalidEventRateMask = ~isfinite(eventRateArray);
    if any(invalidEventRateMask, 'all')
        warning('Session %s contains %d non-finite event-rate values; replacing them with 0.', ...
            outputStem, nnz(invalidEventRateMask));
        eventRateArray(invalidEventRateMask) = 0;
    end

    finiteX = x_position(isfinite(x_position));
    if numel(finiteX) < 2
        error('Session %s does not contain enough finite %s values.', outputStem, opts.xPositionColumn);
    end

    if min(finiteX) == max(finiteX)
        error('Session %s has no %s range for spatial binning.', outputStem, opts.xPositionColumn);
    end

    bin_edges = linspace(min(finiteX), max(finiteX), opts.numBins + 1);
    [counts, ~] = histcounts(x_position, bin_edges);
    if sum(counts) == 0
        error('Session %s produced zero occupancy counts across spatial bins.', outputStem);
    end
    probabilityOfMouseOccupyingBin = counts / sum(counts);

    cellFiringProbabilityPerBin = calculateFiringProbability(eventRateArray, x_position, bin_edges);

    numTopBins = min(5, size(cellFiringProbabilityPerBin, 2));
    topBins = zeros(numNeurons, numTopBins);
    for neuron = 1:numNeurons
        [~, topIndices] = maxk(cellFiringProbabilityPerBin(neuron, :), numTopBins);
        topBins(neuron, :) = topIndices;
    end

    topBinsPath = fullfile(outputDir, [outputStem '_topBins.h5']);
    recreateH5(topBinsPath);
    h5create(topBinsPath, '/topBins', size(topBins), 'Datatype', 'double');
    h5write(topBinsPath, '/topBins', topBins);

    neuronFiringProbability = mean(eventRateArray ~= 0, 1);
    [MI_perCell, MI_perCellperBin] = calculateMutualInformation( ...
        cellFiringProbabilityPerBin, neuronFiringProbability, probabilityOfMouseOccupyingBin);

    validSubset = opts.binsubset(opts.binsubset >= 1 & opts.binsubset <= numel(probabilityOfMouseOccupyingBin));
    MI_perCell_subset = [];
    if ~isempty(validSubset)
        MI_perCell_subset = MI_perCellperBin(:, validSubset) * probabilityOfMouseOccupyingBin(validSubset)';
    end

    actualMiPath = fullfile(outputDir, [outputStem '_MI_per_cell_actual.h5']);
    recreateH5(actualMiPath);
    h5create(actualMiPath, '/MI_perCellActual', size(MI_perCell), 'Datatype', 'double');
    if ~isempty(MI_perCell_subset)
        h5create(actualMiPath, '/MI_perCellSubset', size(MI_perCell_subset), 'Datatype', 'double');
    end
    h5create(actualMiPath, '/MI_perCellperBin', size(MI_perCellperBin), 'Datatype', 'double');
    h5create(actualMiPath, '/binOccupancyProbability', size(probabilityOfMouseOccupyingBin), 'Datatype', 'double');

    h5write(actualMiPath, '/MI_perCellActual', MI_perCell);
    if ~isempty(MI_perCell_subset)
        h5write(actualMiPath, '/MI_perCellSubset', MI_perCell_subset);
    end
    h5write(actualMiPath, '/MI_perCellperBin', MI_perCellperBin);
    h5write(actualMiPath, '/binOccupancyProbability', probabilityOfMouseOccupyingBin);

    numCells = size(cellFiringProbabilityPerBin, 1);
    numBins = size(cellFiringProbabilityPerBin, 2);
    chunkSize = max(1, floor(opts.maxChunkElements / max(1, numCells * numBins)));
    chunkSize = min([chunkSize, opts.numShuffles, opts.maxShufflesPerWrite]);

    miResultsPath = fullfile(outputDir, [outputStem '_MI_results.h5']);
    recreateH5(miResultsPath);
    h5create(miResultsPath, '/MI_perCellAllShuffles', [numCells, opts.numShuffles], ...
        'Datatype', 'double', 'ChunkSize', [numCells, chunkSize]);
    h5create(miResultsPath, '/MI_perCellperBinAllShuffles', [numCells, numBins, opts.numShuffles], ...
        'Datatype', 'double', 'ChunkSize', [numCells, numBins, chunkSize]);

    for startShuffle = 1:chunkSize:opts.numShuffles
        endShuffle = min(startShuffle + chunkSize - 1, opts.numShuffles);
        currentChunkSize = endShuffle - startShuffle + 1;

        MI_perCellperBinAllShufflesChunk = zeros(numCells, numBins, currentChunkSize);
        MI_perCellAllShufflesChunk = zeros(numCells, currentChunkSize);

        for shuffleIdx = 1:currentChunkSize
            shuffle = startShuffle + shuffleIdx - 1;
            disp(['Shuffle ' num2str(shuffle) ' for ' outputStem]);

            shuffledEventRate = shuffleNeuronFrames(eventRateArray, numFrames, numNeurons);
            cellFiringProbabilityPerBinThisShuffle = calculateFiringProbability(shuffledEventRate, x_position, bin_edges);
            neuronFiringProbabilityThisShuffle = mean(shuffledEventRate ~= 0, 1);

            [MI_perCellThisShuffle, MI_perCellperBinThisShuffle] = calculateMutualInformation( ...
                cellFiringProbabilityPerBinThisShuffle, neuronFiringProbabilityThisShuffle, probabilityOfMouseOccupyingBin);

            MI_perCellperBinAllShufflesChunk(:, :, shuffleIdx) = MI_perCellperBinThisShuffle;
            MI_perCellAllShufflesChunk(:, shuffleIdx) = MI_perCellThisShuffle(:);
        end

        h5write(miResultsPath, '/MI_perCellAllShuffles', MI_perCellAllShufflesChunk, ...
            [1, startShuffle], [numCells, currentChunkSize]);
        h5write(miResultsPath, '/MI_perCellperBinAllShuffles', MI_perCellperBinAllShufflesChunk, ...
            [1, 1, startShuffle], [numCells, numBins, currentChunkSize]);
    end
end

function shuffledMatrix = shuffleNeuronFrames(eventRateArray, numFrames, numNeurons)
    shuffledMatrix = eventRateArray;
    for neuron = 1:numNeurons
        shuffleIndices = randperm(numFrames);
        shuffledMatrix(:, neuron) = eventRateArray(shuffleIndices, neuron);
    end
end

function recreateH5(h5FilePath)
    if isfile(h5FilePath)
        delete(h5FilePath);
    end
end

function sessionStem = buildSessionStem(inputBase, sessionIdx, sessionSource, runTimestamp)
    [~, sessionBase, ~] = fileparts(char(sessionSource));
    if isempty(sessionBase)
        sessionBase = ['session_' num2str(sessionIdx)];
    end

    sanitizedSessionBase = regexprep(sessionBase, '[^A-Za-z0-9]+', '_');
    sanitizedSessionBase = regexprep(sanitizedSessionBase, '_+', '_');
    sanitizedSessionBase = regexprep(sanitizedSessionBase, '^_|_$', '');

    if isempty(sanitizedSessionBase)
        sanitizedSessionBase = ['session_' num2str(sessionIdx)];
    end

    maxNameLength = 80;
    if strlength(string(sanitizedSessionBase)) > maxNameLength
        sanitizedSessionBase = extractBefore(string(sanitizedSessionBase), maxNameLength + 1);
        sanitizedSessionBase = char(sanitizedSessionBase);
    end

    sessionStem = sprintf('%s__run%s__session%02d__%s', inputBase, runTimestamp, sessionIdx, sanitizedSessionBase);
end
