%% calc neurons with significant spatial information 
% Define the HDF5 file path and dataset name for shuffled mutual information data
MI_shuffled_results_file = 'dataFinalesLinearTrackMouse6Day4cellTracesAlignedToTracking_MI_results.h5';
h5ResultsFilePath = strcat(filePath, MI_shuffled_results_file);
datasetName = '/MI_perCellAllShuffles';

% Load the shuffled MI data for all 1000 shuffles
MI_perCellAllShuffles = h5read(h5ResultsFilePath, datasetName);

% Compute the 95th percentile (one-sided upper bound) for the shuffled data
upperBound = prctile(MI_perCellAllShuffles, 95, 2); % 95th percentile across shuffles for each cell

%%
% Find indices of neurons where actual MI is significantly greater than the shuffled 95% threshold
significantIndices = find(MI_perCell > upperBound);

% Display the indices of significant neurons
disp('Indices of neurons with significant mutual spatial information:');
disp(significantIndices);

%% get top firing bins of significantIndicies
topBinsSigIndices = topBins(significantIndices, :);

%% for visualiation normalize MI to shuffle
% Compute the mean MI for each cell across all shuffles
meanMI_perCellShuffled = mean(MI_perCellAllShuffles, 2);
% Normalize each cell's actual MI by the mean MI of its shuffled values
normalizedMI_perCell = MI_perCell ./ meanMI_perCellShuffled;
normalizedMI_perCellSignificantIdx = normalizedMI_perCell(significantIndices, :);

%% store 
% Initialize the structure array with the first file's data if needed 
dataStruct(1).fileName = fileName;
dataStruct(1).topBinsSigIndices = topBinsSigIndices;
dataStruct(1).significantIndices = significantIndices;
dataStruct(1).normalizedMI_perCellSignificantIdx = normalizedMI_perCellSignificantIdx;
dataStruct(1).normalizedMI_perCell = normalizedMI_perCell;

%% Suppose you have new values for fileName, topBinsSigIndices, and significantIndices
newIndex = length(dataStruct) + 1;  % Find the next index in the structure array

% Add the new data to the structure array
dataStruct(newIndex).fileName = fileName;
dataStruct(newIndex).topBinsSigIndices = topBinsSigIndices;
dataStruct(newIndex).significantIndices = significantIndices;
dataStruct(newIndex).normalizedMI_perCellSignificantIdx = normalizedMI_perCellSignificantIdx;
dataStruct(newIndex).normalizedMI_perCell = normalizedMI_perCell;

%% load top firing bins 

top_bins_file = 'dataFinalesLinearTrackMouse1Day4cellTracesAlignedToTracking_topBins.h5';
h5topBinsResults = strcat(filePath, top_bins_file);
datasetNameTopBins = '/topBins';

% Load the top bins data 
topBins_perCell = h5read(h5topBinsResults, datasetNameTopBins);

%% vizualize 
% Loop over each neuron and plot the histogram
for neuron = 1:10
    % Get the shuffled distribution for the current neuron
    shuffledData = MI_perCellAllShuffles(neuron, :);
    
    % Create a new figure for each neuron
    figure;
    
    % Plot the histogram of the shuffled data
    histogram(shuffledData, 'Normalization', 'pdf'); % Use 'pdf' to show probability density
    
    % Overlay the actual data value as a vertical line
    hold on;
    xline(MI_perCell(neuron), 'r', 'LineWidth', 2); % Red line for actual value
    hold off;
    
    % Add title and labels
    title(['Neuron ', num2str(neuron), ': Mutual Spatial Information']);
    xlabel('Mutual Spatial Information');
    ylabel('Probability Density');
    legend('Shuffled Distribution', 'Actual Value');
    
    % Optional: pause to view each figure or adjust to save figures as needed
    %pause(0.5); % Pause to allow viewing, or remove if not needed
end
