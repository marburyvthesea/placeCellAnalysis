%% get the number of neurons
cellNumber = size(spikes, 2);

% positionTimes: array of position times (in seconds)
% store time stamps as a single array
positionTimes = pos(1,:);
timeStamps = positionTimes';

% set bin size, etc.
% binSize: size of each spatial bin (e.g., 2.5 cm)
% positions: array of positions corresponding to positionTimes (N x 2 for 2D positions)
pos_transposed = pos';
positions = pos_transposed(:,2:end);
binSize = 3;

% Initialize an array to store the row indices (spike time) for each neuron
spikeTimesArray = cell(1, cellNumber);

% Initialize an array to store the spike positions for each neuron
spikePositionsArray = cell(1, cellNumber);

% Loop through each colume and find the row indices where the value is 1
% spikeTimes: array of spike times (in seconds)
for c = 1:cellNumber
    % Find the row indices where the value in the current colume is 1
    spikeTimes = find(spikes(:, c) == 1);
    % convert row indices to timestamps & store timestamps in the cell array
    spikeTimesArray{c} = timeStamps(spikeTimes);
    spikePositionsArray{c} = positions(spikeTimes, 1);
end


%% function identifyPlaceCells(~, positionTimes, positions, binSize, threshold)

%% Identify place cells by calculating the spatial information score.

    % Define the spatial environment
    % excluding 5% of the total length on both ends of the linear track
    minX = min(positions(:,1)) + 0.05 * (max(positions(:,1)) - min(positions(:,1)));
    maxX = max(positions(:,1)) - 0.05 * (max(positions(:,1)) - min(positions(:,1)));
    minY = min(positions(:,2));
    maxY = max(positions(:,2));

    %% Bin the positions
    edgesX = minX:binSize:maxX;
    %edgesY = minY:binSize:maxY;
    numBinsX = length(edgesX) - 1;
    %numBinsY = length(edgesY) - 1;

    % Calculate occupancy times
    occupancyTimes = zeros(1, numBinsX);
    for i = 1:length(positionTimes)
        x = positions(i, 1);
        %y = positions(i, 2);
        binX = find(x < edgesX, 1) - 1;
        %binY = find(y < edgesY, 1) - 1;
        if ~isempty(binX) && binX > 0
            occupancyTimes(1, binX) = occupancyTimes(1, binX) + 1; % increment time count for the bin
        end
    end

    % Convert occupancy times from counts to seconds
    samplingRate = 1 / mean(diff(positionTimes)); % assuming regular sampling intervals
    % occupancyTimes: a matrix of occupancy times for each spatial bin
    occupancyTimes = occupancyTimes / samplingRate;

%% Calculate spike counts in each bin
% spikeCounts: a matrix of spike counts for each spatial bin
    spikeCounts = zeros(cellNumber, numBinsX);
    for c = 1:cellNumber
        for i = 1:length(spikeTimesArray{c})
            t = spikeTimesArray{c}(i);
            %posIdx = find(positionTimes <= t, 1, 'last'); % find the position closest to the spike time
            posIdx = find(positionTimes == t, 1); % find the position corresponding to the spike time
            x = positions(posIdx, 1);
            binX = find(x < edgesX, 1) - 1;
            if ~isempty(binX) && binX > 0 
                spikeCounts(c, binX) = spikeCounts(c, binX) + 1;
            end
        end
    end

    %% Calculate spatial information score
    % spatialInfoScore = calculateSpatialInformation(spikeCounts, occupancyTimes);

    % Identify place cells
    % if spatialInfoScore > threshold
        % disp(['This neuron is a place cell with spatial information score: ', num2str(spatialInfoScore)]);
    % else
        % disp(['This neuron is not a place cell. Spatial information score: ', num2str(spatialInfoScore)]);
    % end 

%% Function Section: Calculate Spatial Information
% This section contains the function to calculate the spatial information score.
% function spatialInfoScore = calculateSpatialInformation(spikeCounts, occupancyTimes)

%% Calculate the spatial information score

    % Total number of spikes
    totalSpikes = sum(spikeCounts, 2);

    % Total occupancy time
    totalOccupancyTime = sum(occupancyTimes(:));

    % Mean firing rate across all bins
    meanFiringRate = totalSpikes / totalOccupancyTime;

    % Probability of being in each bin
    occupancyProbability = occupancyTimes / totalOccupancyTime;

    % Firing rate in each bin
    firingRate = spikeCounts ./ occupancyTimes;

    % Avoid division by zero
    firingRate(isnan(firingRate)) = 0;

    % Spatial information calculation
    spatialInfoScore = zeros(cellNumber, 1);
    for n = 1:size(spikeCounts, 1)
        for i = 1:size(spikeCounts, 2)
            if occupancyProbability(i) > 0 && firingRate(n, i) > 0
                spatialInfoScore(n, 1) = spatialInfoScore(n, 1) + ...
                    occupancyProbability(i) * (firingRate(n, i) / meanFiringRate(n)) * log2(firingRate(n, i) / meanFiringRate(n));
            end
        end
    end

    % Convert from bits/bin(bits per event?) to bits/spike
    % spatialInfoScore = spatialInfoScore / meanFiringRate;


%% Define place cells
% shuffle event locations for each neuron 1000 times & recalculate spatial
% information for each shuffle

% Number of shuffles
numShuffles = 1000;

% Initialize a cell array to store the shuffled positions
shuffledPositions = cell(cellNumber, numShuffles);

% Loop through each neuron
for n = 1:cellNumber
    % Get the number of events for this neuron
    numEvents = length(spikePositionsArray{n});

    % Perform the shuffling 1000 times
    for shuffleIdx = 1 : numShuffles
        % randomly assign spike positions
        random_indices = randperm(size(positions, 1), numEvents);
        shuffledPositions{n, shuffleIdx} = positions(random_indices, 1);
    end
end


% Calculate spatial information for shuffled event positions
shuffledSpatialInfoScore = zeros(cellNumber, numShuffles);
for n = 1:cellNumber
    % Calculate shuffled spike counts in each bin 
    shuffledSpikeCounts = zeros(numShuffles, numBinsX);
    for s = 1:numShuffles
        for p = 1:length(shuffledPositions{n, s})
            posX = shuffledPositions{n, s}(p);
            binX = find(posX < edgesX, 1) - 1;
            if ~isempty(binX) && binX > 0
                shuffledSpikeCounts(s, binX) = shuffledSpikeCounts(s, binX) + 1;
            end
        end
        shuffledFiringRate = shuffledSpikeCounts ./occupancyTimes;  
        shuffledFiringRate(isnan(shuffledFiringRate)) = 0;
        % Spatial informatioin calculation for shuffled spike positions
        for b = 1:numBinsX
            if occupancyProbability(b) > 0 && shuffledFiringRate(s, b) > 0
                shuffledSpatialInfoScore(n, s) = shuffledSpatialInfoScore(n, s) + ...
                    occupancyProbability(b) * (shuffledFiringRate(s, b) / meanFiringRate(n)) * log2(shuffledFiringRate(s, b) / meanFiringRate(n));
            end
        end
    end
end

% sort each row from smallest number to largest number
shuffledSpatialInfoScore_sorted = sort(shuffledSpatialInfoScore, 2);

% threshold: spatial information score threshold to identify place cells
% get the value at 95% of each shuffled distribution and set as the threshold for place cells
threshold = shuffledSpatialInfoScore_sorted(:, numShuffles * 0.95);

placeCellIdx = find(spatialInfoScore > threshold);


