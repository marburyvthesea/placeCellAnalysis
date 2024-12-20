%% calc neurons with significant spatial information 
% Define the HDF5 file path and dataset name for shuffled mutual information data
filePath = '/Users/johnmarshall/Documents/Analysis/miniscope_lineartrack/mIAnalysis/day2/' ; 
MI_shuffled_results_file = 'dataFinalesLinearTrackMouse3day2cellTracesAlignedToTracking_MI_results.h5';
trackingFile = 'ell_traces_Mouse3day2cellTracesAlignedToTracking.csv';
h5ResultsFilePath = strcat(filePath, MI_shuffled_results_file);
datasetName = '/MI_perCellAllShuffles';

% Load the shuffled MI data for all 1000 shuffles
MI_perCellAllShuffles = h5read(h5ResultsFilePath, datasetName);

% define bins to include 
binsubset = 2:31;
MI_perCellPerBinShuffles = h5read(h5ResultsFilePath, '/MI_perCellperBinAllShuffles');
MI_perCellPerBinShuffles_subset = MI_perCellPerBinShuffles(:, binsubset, :);

% Compute the 95th percentile (one-sided upper bound) for the shuffled data
upperBound = prctile(MI_perCellAllShuffles, 95, 2); % 95th percentile across shuffles for each cell
%%
MI_actual_results_file = 'ell_traces_Mouse3day2_MI_per_cell_actual.h5';
topBins_file = 'day2Mouse3day2cellTracesAlignedToTracking_topBins.h5';
h5ResultsFilePathActual = strcat(filePath, MI_actual_results_file);
MI_perCell = h5read(h5ResultsFilePathActual, '/MI_perCellActual');
MI_perCell_subset = h5read(h5ResultsFilePathActual, '/MI_perCellSubset');
MI_perCellPerBin = h5read(h5ResultsFilePathActual, '/MI_perCellperBin');
binOccupancyProbability = h5read(h5ResultsFilePathActual, '/binOccupancyProbability');

toBinsFilePath = strcat(filePath, topBins_file);
topBins = h5read(toBinsFilePath, '/topBins');

%
MI_perCellShuffled_subset_precat = arrayfun(@(z) MI_perCellPerBinShuffles_subset(:, :, z) * binOccupancyProbability(binsubset)' ...
    , 1:size(MI_perCellPerBinShuffles_subset, 3), 'UniformOutput', false);
MI_perCellShuffled_subset = squeeze(cat(3, MI_perCellShuffled_subset_precat{:}));  % Convert the cell array to a 3D matrix

upperBoundSubset = prctile(MI_perCellShuffled_subset, 95, 2); % 95th percentile across shuffles for each cell
%%
% Find indices of neurons where actual MI is significantly greater than the shuffled 95% threshold
significantIndices = find(MI_perCell > upperBound);
significantIndicesSubset = find(MI_perCell_subset > upperBoundSubset);

% Display the indices of significant neurons
disp('Indices of neurons with significant mutual spatial information:');
disp(significantIndices);

disp('Indices of neurons with significant mutual spatial information subset:');
disp(significantIndicesSubset);

%% get top firing bins of significantIndicies
topBinsSigIndices = topBins(significantIndices, :);

topBinsSigIndicesSubset = topBins(significantIndicesSubset, :);
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
    dataStruct(1).fileName = MI_shuffled_results_file;
    dataStruct(1).trackingFileName = trackingFile;
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
%%
%% store 
i = length(dataStruct)+1; 

dataStruct(i).fileName = MI_shuffled_results_file;
dataStruct(i).trackingFileName = trackingFile;
dataStruct(i).topBinsSigIndices = topBinsSigIndices;
dataStruct(i).significantIndices = significantIndices;
dataStruct(i).normalizedMI_perCellSignificantIdx = normalizedMI_perCellSignificantIdx;
dataStruct(i).normalizedMI_perCell = normalizedMI_perCell;

dataStruct(i).topBinsSigIndicesSubset = topBinsSigIndicesSubset;
dataStruct(i).significantIndicesSubset = significantIndicesSubset;
dataStruct(i).normalizedMI_perCellSignificantIdxSubset = normalizedMI_perCellSignificantIdxSubset;
dataStruct(i).normalizedMI_perCellSubset = normalizedMI_perCellSubset;


%% save 
fileName = strcat(datestr(now, 'mm_dd_yy_HH_MM_SS'), '_MIdata.mat');
matFilePath = fullfile(filePath, fileName);

% Save dataStruct to the specified .mat file
save(matFilePath, 'dataStruct');


