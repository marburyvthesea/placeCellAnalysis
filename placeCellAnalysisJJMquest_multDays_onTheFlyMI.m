%% Split a combined multi-session GCAMP dataframe and run place cell analysis per session.
% This variant avoids writing the large shuffled-peaks cube to disk.
% For each shuffle it computes:
%   shuffled traces -> signalPeaks -> event rate -> shuffled MI
% and writes only the MI outputs.
%
% Expected input:
%   alignedFile = path to a combined CSV like GCAMP_with_velocity.csv
%
% Optional caller overrides:
%   samplingRateHz = miniscope sampling rate in Hz (default: 20)
%   numShuffles = number of shuffle iterations (default: 1000)
%   numBins = number of spatial bins (default: 32)
%   binsubset = subset of bins used for MI subset output (default: 2:31)
%   maxChunkElements = target maximum elements per shuffle-results chunk (default: 5e8)
%   maxShufflesPerWrite = maximum shuffle count to accumulate before writing MI results (default: 100)
%   writeSplitSessionTables = whether to write one CSV per session (default: true)
%   sessionIdColumn = column used to split sessions (default: 'ezTrackOutput')
%   frameColumn = column used to sort rows within session (default: 'closestBehavCamFrameIdx')

if ~exist('alignedFile', 'var') || isempty(alignedFile)
    error('Set alignedFile to the combined GCAMP CSV path before running this script.');
end

if ~exist('samplingRateHz', 'var') || isempty(samplingRateHz)
    samplingRateHz = 20;
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

disp('loading');
disp(alignedFile);
dataTable = readtable(alignedFile, 'VariableNamingRule', 'preserve');

variableNames = dataTable.Properties.VariableNames;
cellColumns = variableNames(startsWith(variableNames, 'cell_'));

if isempty(cellColumns)
    error('No cell_* columns were found in %s.', alignedFile);
end

requiredColumns = {sessionIdColumn, 'X_coor', 'Y_coor'};
for colIdx = 1:numel(requiredColumns)
    if ~ismember(requiredColumns{colIdx}, variableNames)
        error('Required column "%s" was not found in %s.', requiredColumns{colIdx}, alignedFile);
    end
end

sessionIds = string(dataTable.(sessionIdColumn));
validRows = ~ismissing(sessionIds) & strlength(strtrim(sessionIds)) > 0;

if ~all(validRows)
    warning('Dropping %d rows with missing session identifiers.', sum(~validRows));
    dataTable = dataTable(validRows, :);
    sessionIds = sessionIds(validRows);
end

[inputDir, inputBase, ~] = fileparts(alignedFile);
outputDir = fullfile(inputDir, [inputBase '_placeCellAnalysis_multDays_onTheFlyMI']);
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
opts.samplingRateHz = samplingRateHz;
opts.numShuffles = numShuffles;
opts.numBins = numBins;
opts.binsubset = binsubset;
opts.maxChunkElements = maxChunkElements;
opts.maxShufflesPerWrite = maxShufflesPerWrite;
opts.numStdsForThresh = 2.5;

for sessionIdx = 1:numSessions
    sessionSource = sessionSources(sessionIdx);
    sessionMask = sessionIds == sessionSource;
    sessionTable = dataTable(sessionMask, :);

    if ismember(frameColumn, sessionTable.Properties.VariableNames)
        sessionTable = sortrows(sessionTable, frameColumn);
    end

    sessionStem = buildSessionStem(inputBase, sessionIdx, sessionSource);
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

manifestPath = fullfile(outputDir, [inputBase '_session_manifest.csv']);
writetable(manifest, manifestPath);

disp('multi-session place cell analysis completed');
disp(outputDir);

function analyzeSingleSession(sessionTable, cellColumns, outputDir, outputStem, opts)
    cellTracesArray = table2array(sessionTable(:, cellColumns));
    x_position = sessionTable.X_coor;
    y_position = sessionTable.Y_coor;

    [numFrames, numNeurons] = size(cellTracesArray);
    if numFrames < 2
        error('Session %s has fewer than 2 frames.', outputStem);
    end

    if numNeurons < 1
        error('Session %s has no cell columns.', outputStem);
    end

    finiteX = x_position(isfinite(x_position));
    finiteY = y_position(isfinite(y_position));

    if numel(finiteX) < 2
        error('Session %s does not contain enough finite X_coor values.', outputStem);
    end

    if numel(finiteY) < 2
        error('Session %s does not contain enough finite Y_coor values.', outputStem);
    end

    if min(finiteX) == max(finiteX)
        error('Session %s has no X_coor range for spatial binning.', outputStem);
    end

    % Retain a 1-second window on the 20 Hz miniscope-aligned rows used in the notebook.
    window_size = max(1, min(numFrames, round(opts.samplingRateHz)));

    %% calculate 2d velocity
    time = (0:numFrames-1)' ./ opts.samplingRateHz;
    dt = diff(time);
    x_velocity = [NaN; diff(x_position) ./ dt];
    y_velocity = [NaN; diff(y_position) ./ dt];
    velocity_2d = sqrt(x_velocity.^2 + y_velocity.^2); %#ok<NASGU>

    %% calculate spike rate using 1 second sliding window
    signalPeaks = computeSignalPeaks(cellTracesArray', 'doMovAvg', 0, 'reportMidpoint', 1, 'numStdsForThresh', opts.numStdsForThresh)';
    spikes = signalPeaks;
    window = ones(window_size, 1);
    event_rate = conv2(spikes, window, 'same') / window_size;

    %% compute mutual information
    bin_edges = linspace(min(finiteX), max(finiteX), opts.numBins + 1);
    [counts, ~] = histcounts(x_position, bin_edges);
    if sum(counts) == 0
        error('Session %s produced zero occupancy counts across spatial bins.', outputStem);
    end
    probabilityOfMouseOccupyingBin = counts / sum(counts);

    cellFiringProbabilityPerBin = calculateFiringProbability(event_rate, x_position, bin_edges);

    %% top firing bins
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

    %% actual mutual information outputs
    neuronFiringProbability = mean(signalPeaks, 1);
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

    %% calculate shuffled mutual information directly in chunks
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

            shuffledCellTraces = cellTracesArray;
            for neuron = 1:numNeurons
                shuffleIndices = randperm(numFrames);
                shuffledCellTraces(:, neuron) = cellTracesArray(shuffleIndices, neuron);
            end

            signalPeaksThisShuffle = computeSignalPeaks(shuffledCellTraces', 'doMovAvg', 0, ...
                'reportMidpoint', 1, 'numStdsForThresh', opts.numStdsForThresh)';
            event_rateThisShuffle = conv2(signalPeaksThisShuffle, window, 'same') / window_size;
            cellFiringProbabilityPerBinThisShuffle = calculateFiringProbability(event_rateThisShuffle, x_position, bin_edges);
            neuronFiringProbabilityThisShuffle = mean(signalPeaksThisShuffle, 1);

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

function recreateH5(h5FilePath)
    if isfile(h5FilePath)
        delete(h5FilePath);
    end
end

function sessionStem = buildSessionStem(inputBase, sessionIdx, sessionSource)
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

    sessionStem = sprintf('%s__session%02d__%s', inputBase, sessionIdx, sanitizedSessionBase);
end
