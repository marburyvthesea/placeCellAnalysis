%%load csv file with peaks and tracking data
%run("loadSortedSession.m")

%%to run on quest
%import path to .csv file to load dataTable
disp('loading');
disp(alignedFile);
dataTable = readtable(alignedFile, 'VariableNamesLine', 1);
dataTable(1, :) = [];
numFrames=size(dataTable,1);

[filePath, fName, ext] = fileparts(alignedFile);
fileName = strcat(fName, ext); 

%%
frameTimes= dataTable(:,1);
%% converting time row to seconds
timeStrings = frameTimes; % Extract the time strings from the table
timeDurations = duration(extractAfter(timeStrings(:,1).Var1, ' days '), 'Format', 'hh:mm:ss.SSSSSS');
pos = seconds(timeDurations);

cellTraces = dataTable(:, 2:end-7);
X_coor= dataTable.X_coor;
Y_coor= dataTable.Y_coor;


%% check if cells have mean firing rate >0.01Hz during periods of movement

%% calculate 2d velocity 
time = pos(:, 1); % Time in seconds
x_position = X_coor; % X coordinate in cm
y_position = Y_coor; % Y coordinate in cm
% Calculate velocities using central difference formula
dt = diff(time);
x_velocity = [NaN; diff(x_position) ./ dt];
y_velocity = [NaN; diff(y_position) ./ dt];
% Compute overall 2D velocity
velocity_2d = sqrt(x_velocity.^2 + y_velocity.^2);
velocity_2d = [velocity_2d; NaN];

%% calculate spike rate using 1 second sliding window
% label peaks exceeding 2.5 SD threshold
[signalPeaks] = computeSignalPeaks(table2array(cellTraces), 'doMovAvg', 0, 'reportMidpoint', 1, 'numStdsForThresh', 2.5);
%% Extract time information and caclulate event rate 
spikes=signalPeaks;
time = pos;
% Define the 1-second window
window_size = round(1 / mean(diff(time)));
% Create a 1D convolution window to sum across the time dimension
window = ones(window_size, 1);
% Perform sliding window sum along the time dimension using convolution
event_rate = conv2(spikes, window, 'same') / window_size;

%% compute mutual information 
%% track position
% by length: bin track into 4 pixel bins
%bin_edges = 0:4:max(x_position)+4;

% by number: calculate the bin edges for 32 equal-length bins (should be ~4
% cm based on cropping)
numBins = 32;
bin_edges = linspace(min(x_position), max(x_position), numBins + 1);

% Compute histogram counts and bin edges
[counts, edges] = histcounts(x_position, bin_edges);
% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
% Calculate probability of the object occupying each bin
probabilityOfMouseOccupyingBin = counts / sum(counts);


%% firing probability in bins 
% Assuming 'event_rate' is a 98x11839 array and 'x_position' is a 1x11858 array
event_rate_t = event_rate'; 
x_position_t = x_position';

% Initialize a matrix to store the probabilities for each element
cellFiringProbabilityPerBin = zeros(size(event_rate_t, 1), length(counts));

% Loop over each element in 'event_rate'
for i = 1:size(event_rate_t, 1)
    % get indicies where cell event rate is nonzero 
    nonzero_indices = find(event_rate_t(i, :) ~= 0);
    
    %get x location at these indicies 
    x_position_when_active=x_position_t(nonzero_indices);

    % get the counts for when the mouse is active in each bin, in # of
    % samples
    [event_counts_per_x_bin, edges] = histcounts(x_position_when_active, bin_edges);

    % Calculate the probability for each bin
    probability_per_bin = event_counts_per_x_bin ./ length(event_rate_t(i,:));
    cellFiringProbabilityPerBin(i, :) = probability_per_bin;
end

% for each cell in the session, return the bin with highest firing probablities, (top 5) 
numNeurons = size(cellFiringProbabilityPerBin, 1);
topBins = zeros(numNeurons, 5); % Preallocate array to store the top 5 bin indices for each neuron
for neuron = 1:numNeurons
    % Get the top 5 bin indices for the current neuron
    [~, topIndices] = maxk(cellFiringProbabilityPerBin(neuron, :), 5);
    % Store the indices in the result array
    topBins(neuron, :) = topIndices;
end

result = regexp(fileName, '_([^_]+)\.csv$', 'tokens', 'once');
h5FilePath = strcat(filePath,'/', fName(1:15),result{1}, '_topBins.h5');  % Set the HDF5 file path
datasetName = '/topBins';  % Define the dataset name
if isfile(h5FilePath)
    delete(h5FilePath);
end
h5create(h5FilePath, datasetName, size(topBins), 'Datatype', 'double');
h5write(h5FilePath, datasetName, topBins);

%% calculate mutual information per cell

% for subset calulation
binsubset = 2:31;

% overall firing probability for each neuron 
neuronFiringProbability = mean(signalPeaks, 1);

% probability of calcium event over occupany probability
expanded_probability_per_bin = repmat(probabilityOfMouseOccupyingBin, size(cellFiringProbabilityPerBin, 1), 1);
p_event_over_p_occupancy = cellFiringProbabilityPerBin ./ expanded_probability_per_bin;

% sum of event probabilities in each bin
MI_perCellperBin = zeros(size(cellFiringProbabilityPerBin,1), size(cellFiringProbabilityPerBin,2));
MI_perCell = zeros(1, size(cellFiringProbabilityPerBin,1));

% for each cell
size(cellFiringProbabilityPerBin,1)
for i = 1:size(cellFiringProbabilityPerBin,1)

    MI_perbin = zeros(1, size(probabilityOfMouseOccupyingBin, 2));
    
    %for each bin  
    for j = 1:size(probabilityOfMouseOccupyingBin, 2)

        %%TO DO: from here 

        % conditional p of event given in bin is = (prob in bin + prob of event) / (prob in bin)

        % (conditional probability of observing a calcium events given mouse is in bin) *    
        % log ( (conditional probability of observing a calcium events given mouse is in bin) / (probability of observing k calcium events (0 or 1) for a cell
        jointProbability = neuronFiringProbability(1, i) .* probabilityOfMouseOccupyingBin(1, j);
        cProbEventGivenBin = jointProbability/probabilityOfMouseOccupyingBin(1, j); 
        
        SI_thisbin = cProbEventGivenBin*log(cProbEventGivenBin/neuronFiringProbability(1, i));
        
        MI_perbin(1, j) = SI_thisbin; 

    end

    MI_perCellperBin(i,:) = MI_perbin; 

MI_perCell = MI_perCellperBin * probabilityOfMouseOccupyingBin'; 

MI_perCell_subset = MI_perCellperBin(:, binsubset) * probabilityOfMouseOccupyingBin(binsubset)';
     
end

h5FilePath = strcat(filePath,'/', fName(1:21), '_MI_per_cell_actual.h5');  % Set the HDF5 file path

if isfile(h5FilePath)
    delete(h5FilePath);
end
h5create(h5FilePath, '/MI_perCellActual', size(MI_perCell), 'Datatype', 'double');
h5create(h5FilePath, '/MI_perCellSubset', size(MI_perCell_subset), 'Datatype', 'double');
h5create(h5FilePath, '/MI_perCellperBin', size(MI_perCellperBin), 'Datatype', 'double');
h5create(h5FilePath, '/binOccupancyProbability', size(probabilityOfMouseOccupyingBin), 'Datatype', 'double');

h5write(h5FilePath, '/MI_perCellActual', MI_perCell);
h5write(h5FilePath, '/MI_perCellSubset', MI_perCell_subset);
h5write(h5FilePath, '/MI_perCellperBin', MI_perCellperBin);
h5write(h5FilePath, '/binOccupancyProbability', probabilityOfMouseOccupyingBin);

%% shuffle calcium traces and compute signal peaks in chunks to avoid filling up RAM

% Define parameters
numShuffles = 1000;
chunkSize = floor(500000000 / (numFrames*numNeurons)); % Number of shuffles per chunk
cellTracesArray = table2array(cellTraces);
[numFrames, numNeurons] = size(cellTracesArray);
% Define HDF5 file path and dataset name
result = regexp(fileName, '_([^_]+)\.csv$', 'tokens', 'once');
h5FilePath = strcat(filePath, '/', fName(1:21), result{1}, '_allPeaksShuffled.h5');
datasetName = '/allPeaksShuffled';


% Initialize the HDF5 dataset with unlimited third dimension
if isfile(h5FilePath)
    delete(h5FilePath); % Remove file if it already exists
end

h5create(h5FilePath, datasetName, [numFrames, numNeurons, Inf], 'Datatype', 'double', 'ChunkSize', [numFrames, numNeurons, chunkSize]);

% Loop over shuffles in chunks
for startShuffle = 1:chunkSize:numShuffles
    % Determine the end of the current chunk
    endShuffle = min(startShuffle + chunkSize - 1, numShuffles);
    currentChunkSize = endShuffle - startShuffle + 1;

    % Initialize array for the current chunk of shuffles
    allPeaksShuffledChunk = zeros(numFrames, numNeurons, currentChunkSize);

    % Perform shuffles in the current chunk
    for shuffle = 1:currentChunkSize
        disp(['Shuffle ' num2str(startShuffle + shuffle - 1)]);

        % Initialize an array to store the shuffled data
        shuffledCellTraces = cellTracesArray;

        % Loop through each neuron (column) and shuffle the frames
        for neuron = 1:numNeurons
            % Generate a random permutation of row indices
            shuffleIndices = randperm(numFrames);
            % Apply the shuffle to the current neuron (column)
            shuffledCellTraces(:, neuron) = cellTracesArray(shuffleIndices, neuron);
        end

        % Compute signal peaks for this shuffle
        signalPeaksThisShuffle = computeSignalPeaks(shuffledCellTraces, 'doMovAvg', 0, 'reportMidpoint', 1, 'numStdsForThresh', 2.5);
        
        % Store the result in the current chunk
        allPeaksShuffledChunk(:, :, shuffle) = signalPeaksThisShuffle;
    end
    % Define the count of elements to write, matching the size of allPeaksShuffledChunk
    count = [numFrames, numNeurons, currentChunkSize];
    % Write the current chunk to the HDF5 file
    h5write(h5FilePath, datasetName, allPeaksShuffledChunk, [1, 1, startShuffle], count);
end

disp('All shuffles have been processed and saved to the HDF5 file.');

%% calculate mutual information iteratively on chunks

% Define HDF5 file paths
h5FilePath = strcat(filePath,'/', fName(1:21), result{1}, '_allPeaksShuffled.h5');  % Path for shuffled data
h5ResultsFilePath = strcat(filePath,'/', fName(1:21),result{1}, '_MI_results.h5');  % Path for MI results

% Define the size of each dataset
[numCells, numBins] = size(cellFiringProbabilityPerBin);  % Assuming dimensions of cells and bins
numShuffles = 1000;  % Total number of shuffles
chunkSize = 100;  % Number of shuffles per chunk

% Initialize HDF5 datasets for storing mutual information results
% MI_perCellAllShuffles: Stores mutual information for each cell across shuffles
h5create(h5ResultsFilePath, '/MI_perCellAllShuffles', [numCells, numShuffles], 'Datatype', 'double', 'ChunkSize', [numCells, chunkSize]);

% MI_perCellperBinAllShuffles: Stores mutual information per cell per bin across shuffles
h5create(h5ResultsFilePath, '/MI_perCellperBinAllShuffles', [numCells, numBins, numShuffles], 'Datatype', 'double', 'ChunkSize', [numCells, numBins, chunkSize]);

for startShuffle = 1:chunkSize:numShuffles
    % Define the end of the current chunk
    endShuffle = min(startShuffle + chunkSize - 1, numShuffles);
    currentChunkSize = endShuffle - startShuffle + 1;

    % Load shuffled data from the HDF5 file for the current chunk
    allPeaksShuffledChunk = h5read(h5FilePath, '/allPeaksShuffled', [1, 1, startShuffle], [size(cellTracesArray, 1), size(cellTracesArray, 2), currentChunkSize]);

    % Initialize arrays to store mutual information results for the current chunk
    MI_perCellperBinAllShufflesChunk = zeros(numCells, numBins, currentChunkSize);
    MI_perCellAllShufflesChunk = zeros(numCells, currentChunkSize);

    % Loop through each shuffle in the current chunk
    for shuffleIdx = 1:currentChunkSize
        shuffle = startShuffle + shuffleIdx - 1;
        disp(['Calculating mutual information for shuffle: ', num2str(shuffle)]);

        % Extract peaks for the current shuffle
        spikes = allPeaksShuffledChunk(:, :, shuffleIdx);

        % Perform mutual information calculations (as in your code)
        % Calculate event rate
        disp('Calculating event rate');
        event_rateThisShuffle = conv2(spikes, ones(round(1 / mean(diff(pos))), 1), 'same') / round(1 / mean(diff(pos)));

        % Calculate firing probability in bins
        cellFiringProbabilityPerBinThisShuffle = calculateFiringProbability(event_rateThisShuffle, x_position, bin_edges);

        % Calculate overall firing probability per neuron
        neuronFiringProbabilityThisShuffle = mean(spikes, 1);

        % Calculate mutual information for each neuron across bins
        disp('Calculating mutual information per cell');
        [MI_perCellThisShuffle, MI_perCellperBinThisShuffle] = calculateMutualInformation(cellFiringProbabilityPerBinThisShuffle, neuronFiringProbabilityThisShuffle, probabilityOfMouseOccupyingBin);

        % Store results for this shuffle in the current chunk
        MI_perCellperBinAllShufflesChunk(:, :, shuffleIdx) = MI_perCellperBinThisShuffle;
        MI_perCellAllShufflesChunk(:, shuffleIdx) = MI_perCellThisShuffle;
    end

    % Write the current chunk of results to the HDF5 file
    h5write(h5ResultsFilePath, '/MI_perCellAllShuffles', MI_perCellAllShufflesChunk, [1, startShuffle], [numCells, currentChunkSize]);
    h5write(h5ResultsFilePath, '/MI_perCellperBinAllShuffles', MI_perCellperBinAllShufflesChunk, [1, 1, startShuffle], [numCells, numBins, currentChunkSize]);
end

disp('Mutual information calculations and storage completed.');
















