function placeCellAnalysisJJMquest_multDays_savedShuffles(alignedFile, config)
% Shared multi-session place-cell workflow for saved shuffled-peak variants.
% Supported storage modes:
%   'compressed_uint8'    - full shuffled peaks saved as compressed uint8 HDF5
%   'sparse_linear_idx'   - shuffled peaks saved as sparse linear indices in HDF5

    if nargin < 1 || isempty(alignedFile)
        error('Set alignedFile to the combined GCAMP CSV path before running this function.');
    end

    if nargin < 2
        config = struct();
    end

    opts = applyDefaults(config);

    disp('loading');
    disp(alignedFile);
    dataTable = readtable(alignedFile, 'VariableNamingRule', 'preserve');

    variableNames = dataTable.Properties.VariableNames;
    cellColumns = variableNames(startsWith(variableNames, 'cell_'));

    if isempty(cellColumns)
        error('No cell_* columns were found in %s.', alignedFile);
    end

    requiredColumns = {opts.sessionIdColumn, 'X_coor', 'Y_coor'};
    for colIdx = 1:numel(requiredColumns)
        if ~ismember(requiredColumns{colIdx}, variableNames)
            error('Required column "%s" was not found in %s.', requiredColumns{colIdx}, alignedFile);
        end
    end

    sessionIds = string(dataTable.(opts.sessionIdColumn));
    validRows = ~ismissing(sessionIds) & strlength(strtrim(sessionIds)) > 0;

    if ~all(validRows)
        warning('Dropping %d rows with missing session identifiers.', sum(~validRows));
        dataTable = dataTable(validRows, :);
        sessionIds = sessionIds(validRows);
    end

    [inputDir, inputBase, ~] = fileparts(alignedFile);
    outputDir = fullfile(inputDir, [inputBase opts.outputSuffix '_' opts.runTimestamp]);
    splitDir = fullfile(outputDir, 'splitSessions');

    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    if opts.writeSplitSessionTables && ~isfolder(splitDir)
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

    for sessionIdx = 1:numSessions
        sessionSource = sessionSources(sessionIdx);
        sessionMask = sessionIds == sessionSource;
        sessionTable = dataTable(sessionMask, :);

        if ismember(opts.frameColumn, sessionTable.Properties.VariableNames)
            sessionTable = sortrows(sessionTable, opts.frameColumn);
        end

        sessionStem = buildSessionStem(inputBase, sessionIdx, sessionSource, opts.runTimestamp);
        splitCsvPath = "";

        manifest.sessionIndex(sessionIdx) = sessionIdx;
        manifest.sessionSource(sessionIdx) = sessionSource;
        manifest.numRows(sessionIdx) = height(sessionTable);
        manifest.sessionStem(sessionIdx) = sessionStem;

        disp(['processing session ' num2str(sessionIdx) ' of ' num2str(numSessions)]);
        disp(sessionSource);

        try
            if opts.writeSplitSessionTables
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

    manifestPath = fullfile(outputDir, [inputBase '__run' opts.runTimestamp '_session_manifest.csv']);
    writetable(manifest, manifestPath);

    disp('multi-session place cell analysis completed');
    disp(outputDir);
end

function opts = applyDefaults(config)
    opts = config;

    if ~isfield(opts, 'samplingRateHz') || isempty(opts.samplingRateHz)
        opts.samplingRateHz = 20;
    end
    if ~isfield(opts, 'numShuffles') || isempty(opts.numShuffles)
        opts.numShuffles = 1000;
    end
    if ~isfield(opts, 'numBins') || isempty(opts.numBins)
        opts.numBins = 32;
    end
    if ~isfield(opts, 'binsubset') || isempty(opts.binsubset)
        opts.binsubset = 2:31;
    end
    if ~isfield(opts, 'maxChunkElements') || isempty(opts.maxChunkElements)
        opts.maxChunkElements = 500000000;
    end
    if ~isfield(opts, 'writeSplitSessionTables') || isempty(opts.writeSplitSessionTables)
        opts.writeSplitSessionTables = true;
    end
    if ~isfield(opts, 'sessionIdColumn') || isempty(opts.sessionIdColumn)
        opts.sessionIdColumn = 'ezTrackOutput';
    end
    if ~isfield(opts, 'frameColumn') || isempty(opts.frameColumn)
        opts.frameColumn = 'closestBehavCamFrameIdx';
    end
    if ~isfield(opts, 'numStdsForThresh') || isempty(opts.numStdsForThresh)
        opts.numStdsForThresh = 2.5;
    end
    if ~isfield(opts, 'peakDeflateLevel') || isempty(opts.peakDeflateLevel)
        opts.peakDeflateLevel = 9;
    end
    if ~isfield(opts, 'peakChunkShuffles') || isempty(opts.peakChunkShuffles)
        opts.peakChunkShuffles = 1;
    end
    if ~isfield(opts, 'indexChunkEntries') || isempty(opts.indexChunkEntries)
        opts.indexChunkEntries = 1000000;
    end
    if ~isfield(opts, 'storageMode') || isempty(opts.storageMode)
        error('config.storageMode is required.');
    end
    if ~isfield(opts, 'outputSuffix') || isempty(opts.outputSuffix)
        error('config.outputSuffix is required.');
    end
    if ~isfield(opts, 'runTimestamp') || isempty(opts.runTimestamp)
        opts.runTimestamp = datestr(now, 'yyyymmdd_HHMMSS');
    end
end

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

    time = (0:numFrames-1)' ./ opts.samplingRateHz;
    dt = diff(time);
    x_velocity = [NaN; diff(x_position) ./ dt];
    y_velocity = [NaN; diff(y_position) ./ dt];
    velocity_2d = sqrt(x_velocity.^2 + y_velocity.^2); %#ok<NASGU>

    window_size = max(1, min(numFrames, round(opts.samplingRateHz)));
    window = ones(window_size, 1);

    signalPeaks = computeSignalPeaks(cellTracesArray', 'doMovAvg', 0, ...
        'reportMidpoint', 1, 'numStdsForThresh', opts.numStdsForThresh)';
    spikes = signalPeaks;
    event_rate = conv2(spikes, window, 'same') / window_size;

    bin_edges = linspace(min(finiteX), max(finiteX), opts.numBins + 1);
    [counts, ~] = histcounts(x_position, bin_edges);
    if sum(counts) == 0
        error('Session %s produced zero occupancy counts across spatial bins.', outputStem);
    end
    probabilityOfMouseOccupyingBin = counts / sum(counts);
    cellFiringProbabilityPerBin = calculateFiringProbability(event_rate, x_position, bin_edges);

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

    shuffleWriteBatch = max(1, floor(opts.maxChunkElements / max(1, numFrames * numNeurons)));
    shuffleWriteBatch = min(shuffleWriteBatch, opts.numShuffles);

    switch opts.storageMode
        case 'compressed_uint8'
            peaksPath = fullfile(outputDir, [outputStem '_allPeaksShuffled_uint8.h5']);
            writeCompressedPeaks(peaksPath, cellTracesArray, numFrames, numNeurons, outputStem, opts, shuffleWriteBatch);
            computeMiFromCompressedPeaks(peaksPath, outputDir, outputStem, x_position, bin_edges, ...
                probabilityOfMouseOccupyingBin, numFrames, numNeurons, window, window_size, cellFiringProbabilityPerBin, opts);
        case 'sparse_linear_idx'
            peaksPath = fullfile(outputDir, [outputStem '_allPeaksShuffled_sparseIdx.h5']);
            writeSparseLinearIndexPeaks(peaksPath, cellTracesArray, numFrames, numNeurons, outputStem, opts);
            computeMiFromSparseLinearIndexPeaks(peaksPath, outputDir, outputStem, x_position, bin_edges, ...
                probabilityOfMouseOccupyingBin, numFrames, numNeurons, window, window_size, cellFiringProbabilityPerBin, opts);
        otherwise
            error('Unsupported storage mode: %s', opts.storageMode);
    end
end

function writeCompressedPeaks(peaksPath, cellTracesArray, numFrames, numNeurons, outputStem, opts, shuffleWriteBatch)
    recreateH5(peaksPath);
    storageChunkShuffles = min(opts.peakChunkShuffles, opts.numShuffles);
    h5create(peaksPath, '/allPeaksShuffled', [numFrames, numNeurons, Inf], ...
        'Datatype', 'uint8', ...
        'ChunkSize', [numFrames, numNeurons, storageChunkShuffles], ...
        'Deflate', opts.peakDeflateLevel);

    for startShuffle = 1:shuffleWriteBatch:opts.numShuffles
        endShuffle = min(startShuffle + shuffleWriteBatch - 1, opts.numShuffles);
        currentChunkSize = endShuffle - startShuffle + 1;
        allPeaksShuffledChunk = zeros(numFrames, numNeurons, currentChunkSize, 'uint8');

        for shuffle = 1:currentChunkSize
            shuffleNumber = startShuffle + shuffle - 1;
            disp(['Shuffle ' num2str(shuffleNumber) ' for ' outputStem]);

            shuffledCellTraces = shuffleNeuronFrames(cellTracesArray, numFrames, numNeurons);
            signalPeaksThisShuffle = computeSignalPeaks(shuffledCellTraces', 'doMovAvg', 0, ...
                'reportMidpoint', 1, 'numStdsForThresh', opts.numStdsForThresh)';
            allPeaksShuffledChunk(:, :, shuffle) = uint8(signalPeaksThisShuffle);
        end

        count = [numFrames, numNeurons, currentChunkSize];
        h5write(peaksPath, '/allPeaksShuffled', allPeaksShuffledChunk, [1, 1, startShuffle], count);
    end
end

function computeMiFromCompressedPeaks(peaksPath, outputDir, outputStem, x_position, bin_edges, probabilityOfMouseOccupyingBin, numFrames, numNeurons, window, window_size, cellFiringProbabilityPerBin, opts)
    numCells = size(cellFiringProbabilityPerBin, 1);
    numBins = size(cellFiringProbabilityPerBin, 2);
    chunkSize = max(1, floor(opts.maxChunkElements / max(1, numFrames * numNeurons)));
    chunkSize = min(chunkSize, opts.numShuffles);

    miResultsPath = fullfile(outputDir, [outputStem '_MI_results.h5']);
    recreateH5(miResultsPath);
    h5create(miResultsPath, '/MI_perCellAllShuffles', [numCells, opts.numShuffles], ...
        'Datatype', 'double', 'ChunkSize', [numCells, chunkSize]);
    h5create(miResultsPath, '/MI_perCellperBinAllShuffles', [numCells, numBins, opts.numShuffles], ...
        'Datatype', 'double', 'ChunkSize', [numCells, numBins, chunkSize]);

    for startShuffle = 1:chunkSize:opts.numShuffles
        endShuffle = min(startShuffle + chunkSize - 1, opts.numShuffles);
        currentChunkSize = endShuffle - startShuffle + 1;

        allPeaksShuffledChunk = h5read(peaksPath, '/allPeaksShuffled', ...
            [1, 1, startShuffle], [numFrames, numNeurons, currentChunkSize]);

        MI_perCellperBinAllShufflesChunk = zeros(numCells, numBins, currentChunkSize);
        MI_perCellAllShufflesChunk = zeros(numCells, currentChunkSize);

        for shuffleIdx = 1:currentChunkSize
            shuffleNumber = startShuffle + shuffleIdx - 1;
            disp(['Calculating mutual information for shuffle: ' num2str(shuffleNumber) ' for ' outputStem]);

            spikes = double(allPeaksShuffledChunk(:, :, shuffleIdx));
            event_rateThisShuffle = conv2(spikes, window, 'same') / window_size;
            cellFiringProbabilityPerBinThisShuffle = calculateFiringProbability(event_rateThisShuffle, x_position, bin_edges);
            neuronFiringProbabilityThisShuffle = mean(spikes, 1);

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

function writeSparseLinearIndexPeaks(peaksPath, cellTracesArray, numFrames, numNeurons, outputStem, opts)
    recreateH5(peaksPath);

    maxLinearIndex = double(numFrames) * double(numNeurons);
    if maxLinearIndex <= double(intmax('uint32'))
        indexDatatype = 'uint32';
    else
        indexDatatype = 'uint64';
    end

    h5create(peaksPath, '/shuffleStartIdx', [opts.numShuffles, 1], 'Datatype', 'uint64');
    h5create(peaksPath, '/shuffleCounts', [opts.numShuffles, 1], 'Datatype', 'uint64');
    h5create(peaksPath, '/matrixShape', [2, 1], 'Datatype', 'uint64');
    h5create(peaksPath, '/peakLinearIdx', [Inf, 1], ...
        'Datatype', indexDatatype, ...
        'ChunkSize', [opts.indexChunkEntries, 1], ...
        'Deflate', opts.peakDeflateLevel);

    h5write(peaksPath, '/matrixShape', uint64([numFrames; numNeurons]));

    shuffleStartIdx = zeros(opts.numShuffles, 1, 'uint64');
    shuffleCounts = zeros(opts.numShuffles, 1, 'uint64');
    currentOffset = uint64(0);

    for shuffle = 1:opts.numShuffles
        disp(['Shuffle ' num2str(shuffle) ' for ' outputStem]);

        shuffledCellTraces = shuffleNeuronFrames(cellTracesArray, numFrames, numNeurons);
        signalPeaksThisShuffle = computeSignalPeaks(shuffledCellTraces', 'doMovAvg', 0, ...
            'reportMidpoint', 1, 'numStdsForThresh', opts.numStdsForThresh)';

        linearIdx = find(signalPeaksThisShuffle);
        count = uint64(numel(linearIdx));
        shuffleCounts(shuffle) = count;

        if count > 0
            shuffleStartIdx(shuffle) = currentOffset + 1;
            linearIdx = cast(linearIdx(:), indexDatatype);
            h5write(peaksPath, '/peakLinearIdx', linearIdx, [double(currentOffset) + 1, 1], [double(count), 1]);
            currentOffset = currentOffset + count;
        end
    end

    h5write(peaksPath, '/shuffleStartIdx', shuffleStartIdx);
    h5write(peaksPath, '/shuffleCounts', shuffleCounts);
end

function computeMiFromSparseLinearIndexPeaks(peaksPath, outputDir, outputStem, x_position, bin_edges, probabilityOfMouseOccupyingBin, numFrames, numNeurons, window, window_size, cellFiringProbabilityPerBin, opts)
    numCells = size(cellFiringProbabilityPerBin, 1);
    numBins = size(cellFiringProbabilityPerBin, 2);
    chunkSize = max(1, floor(opts.maxChunkElements / max(1, numCells * numBins)));
    chunkSize = min(chunkSize, opts.numShuffles);

    shuffleStartIdx = uint64(h5read(peaksPath, '/shuffleStartIdx'));
    shuffleCounts = uint64(h5read(peaksPath, '/shuffleCounts'));

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
            shuffleNumber = startShuffle + shuffleIdx - 1;
            disp(['Calculating mutual information for shuffle: ' num2str(shuffleNumber) ' for ' outputStem]);

            spikes = zeros(numFrames, numNeurons);
            count = double(shuffleCounts(shuffleNumber));
            if count > 0
                startIdx = double(shuffleStartIdx(shuffleNumber));
                linearIdx = h5read(peaksPath, '/peakLinearIdx', [startIdx, 1], [count, 1]);
                spikes(double(linearIdx)) = 1;
            end

            event_rateThisShuffle = conv2(spikes, window, 'same') / window_size;
            cellFiringProbabilityPerBinThisShuffle = calculateFiringProbability(event_rateThisShuffle, x_position, bin_edges);
            neuronFiringProbabilityThisShuffle = mean(spikes, 1);

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

function shuffledCellTraces = shuffleNeuronFrames(cellTracesArray, numFrames, numNeurons)
    shuffledCellTraces = cellTracesArray;
    for neuron = 1:numNeurons
        shuffleIndices = randperm(numFrames);
        shuffledCellTraces(:, neuron) = cellTracesArray(shuffleIndices, neuron);
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
