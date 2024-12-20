%% calc neurons with significant spatial information 
% Define the HDF5 file path and dataset name for shuffled mutual information data
filePath = '/Users/johnmarshall/Documents/Analysis/miniscope_lineartrack/mIAnalysis/' ; 
MI_shuffled_results_file = 'ell_traces_Mouse1Day2Mouse1Day2cellTracesAlignedToTracking_MI_results.h5';
h5ResultsFilePath = strcat(filePath, MI_shuffled_results_file);
datasetName = '/MI_perCellAllShuffles';

% Load the shuffled MI data for all 1000 shuffles
MI_perCellAllShuffles = h5read(h5ResultsFilePath, datasetName);

% define bins to include 
binsubset = 2:31;
MI_perCellAllShuffles_subset = MI_perCellAllShuffles(:, binsubset, :);

% Compute the 95th percentile (one-sided upper bound) for the shuffled data
upperBound = prctile(MI_perCellAllShuffles, 95, 2); % 95th percentile across shuffles for each cell
upperBoundSubset = prctile(MI_perCellAllShuffles_subset, 95, 2); % 95th percentile across shuffles for each cell
%%
MI_actual_results_file = 'ell_traces_Mouse1Day2_MI_per_cell_actual.h5';
topBins_file = 'ell_traces_MousMouse1Day2cellTracesAlignedToTracking_topBins.h5';
h5ResultsFilePathActual = strcat(filePath, MI_actual_results_file);
MI_perCell = h5read(h5ResultsFilePathActual, '/MI_perCellActual');
MI_perCellperBin = h5read(h5ResultsFilePathActual, '/MI_perCellperBin');

toBinsFilePath = strcat(filePath, topBins_file);
topBins = h5read(toBinsFilePath, '/topBins');

%%
% Find indices of neurons where actual MI is significantly greater than the shuffled 95% threshold
significantIndices = find(MI_perCell > upperBound);
%significantIndicesSubset = find(MI_perCell_subset > upperBoundSubset);

% Display the indices of significant neurons
disp('Indices of neurons with significant mutual spatial information:');
disp(significantIndices);

disp('Indices of neurons with significant mutual spatial information subset:');
disp(significantIndicesSubset);

%% get top firing bins of significantIndicies
topBinsSigIndices = topBins(significantIndices, :);

%topBinsSigIndicesSubset = topBins(significantIndicesSubset, :);
%topBinsSigIndicesSubset(~ismember(topBinsSigIndicesSubset, binsubset)) = NaN;

%% for visualiation normalize MI to shuffle
% Compute the mean MI for each cell across all shuffles
meanMI_perCellShuffled = mean(MI_perCellAllShuffles, 2);
% Normalize each cell's actual MI by the mean MI of its shuffled values
normalizedMI_perCell = MI_perCell ./ meanMI_perCellShuffled;
normalizedMI_perCellSignificantIdx = normalizedMI_perCell(significantIndices, :);

meanMI_perCellShuffledSubset = mean(MI_perCellAllShuffles_subset, 2);
% Normalize each cell's actual MI by the mean MI of its shuffled values
normalizedMI_perCellSubset = MI_perCell_subset ./ meanMI_perCellShuffledSubset;
normalizedMI_perCellSignificantIdxSubset = normalizedMI_perCellSubset(significantIndices, :);

%% store 
% Initialize the structure array with the first file's data if needed 
% Check if dataStruct already exists and is non-empty
if ~exist('dataStruct', 'var') || isempty(dataStruct)
    dataStruct(1).fileName = fileName;
    dataStruct(1).topBinsSigIndices = topBinsSigIndices;
    dataStruct(1).significantIndices = significantIndices;
    dataStruct(1).normalizedMI_perCellSignificantIdx = normalizedMI_perCellSignificantIdx;
    dataStruct(1).normalizedMI_perCell = normalizedMI_perCell;

    dataStruct(1).topBinsSigIndicesSubset = topBinsSigIndicesSubset;
    dataStruct(1).significantIndicesSubset = significantIndicesSubset;
    dataStruct(1).normalizedMI_perCellSignificantIdxSubset = normalizedMI_perCellSignificantIdxSubset;
    dataStruct(1).normalizedMI_perCellSubset = normalizedMI_perCellSubset;

else
    disp('dataStruct already exists; skipping initialization.');
end
%% Suppose you have new values for fileName, topBinsSigIndices, and significantIndices
newIndex = length(dataStruct) + 1;  % Find the next index in the structure array

% Add the new data to the structure array
dataStruct(newIndex).fileName = fileName;
dataStruct(newIndex).topBinsSigIndices = topBinsSigIndices;
dataStruct(newIndex).significantIndices = significantIndices;
dataStruct(newIndex).normalizedMI_perCellSignificantIdx = normalizedMI_perCellSignificantIdx;
dataStruct(newIndex).normalizedMI_perCell = normalizedMI_perCell;

dataStruct(newIndex).topBinsSigIndicesSubset = topBinsSigIndicesSubset;
dataStruct(newIndex).significantIndicesSubset = significantIndicesSubset;
dataStruct(newIndex).normalizedMI_perCellSignificantIdxSubset = normalizedMI_perCellSignificantIdxSubset;
dataStruct(newIndex).normalizedMI_perCellSubset = normalizedMI_perCellSubset;
%% save 
fileName = strcat(datestr(now, 'mm_dd_yy_HH_MM_SS'), '_MIdata.mat');
matFilePath = fullfile(filePath, fileName);

% Save dataStruct to the specified .mat file
save(matFilePath, 'dataStruct');

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
