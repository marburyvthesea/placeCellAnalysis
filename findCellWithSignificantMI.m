%% calculate neurons with significant spatial information

directoryPath = '/Users/johnmarshall/Documents/Analysis/miniscope_analysis/miniscopeLinearTrack/mIAnalysis/';
MI_resultsFileShuffled = '388_Run4_20250707_motcorrectedcellTracesAlignedToTracking_MI_results.h5'; 
MI_resultsFileActual = '388_Run4_20250707_mot_MI_per_cell_actual.h5';
trackingFile = '388_Run4_20250707_motion_correctedcellTracesAlignedToTracking.csv';

h5ResultsFilePathShuffled = strcat(directoryPath, MI_resultsFileShuffled);
h5ResultsFilePathActual = strcat(directoryPath, MI_resultsFileActual);

% full dataset - including ends of track 
shuffledDatasetName = '/MI_perCellAllShuffles';
actualDatasetName = '/MI_perCellActual';
MI_perCellAllShuffles = h5read(h5ResultsFilePathShuffled, shuffledDatasetName);
MI_perCellActual = h5read(h5ResultsFilePathActual, actualDatasetName);

% subset dataset - exclude ends of track 
% for shuffled data use "MI per bin" and exclude bins at end of track  
binsubset = 2:31;
MI_perCellPerBinShuffles = h5read(h5ResultsFilePathShuffled, '/MI_perCellperBinAllShuffles');
MI_perCellPerBinShuffles_subset = MI_perCellPerBinShuffles(:, binsubset, :);

% MI calculation on actual data excluding ends of track 
actualDatasetNameSubset = '/MI_perCellSubset';
MI_perCellActual_subset = h5read(h5ResultsFilePathActual, actualDatasetNameSubset);


%% Compute the 95th percentile (one-sided upper bound) for the shuffled data
% all bins
upperBoundAllShuffles = prctile(MI_perCellAllShuffles, 95, 2); % 95th percentile across shuffles for each cell

% bin subsets
MI_perCellPerBin = h5read(h5ResultsFilePathActual, '/MI_perCellperBin');
binOccupancyProbability = h5read(h5ResultsFilePathActual, '/binOccupancyProbability');
%
MI_perCellShuffled_subset_precat = arrayfun(@(z) MI_perCellPerBinShuffles_subset(:, :, z) * binOccupancyProbability(binsubset)' ...
    , 1:size(MI_perCellPerBinShuffles_subset, 3), 'UniformOutput', false);
MI_perCellShuffled_subset = squeeze(cat(3, MI_perCellShuffled_subset_precat{:}));  % Convert the cell array to a 3D matrix

upperBoundSubset = prctile(MI_perCellShuffled_subset, 95, 2); 

%% Find indices of neurons where actual MI is significantly greater than the shuffled 95% threshold
significantIndices = find(MI_perCellActual > upperBoundAllShuffles);
significantIndicesSubset = find(MI_perCellActual_subset > upperBoundSubset);

% Display the indices of significant neurons - using all bins
disp('Indices of neurons with significant mutual spatial information including track ends:');
disp(significantIndices);

% Display the indices of significant neurons - excluding track ends 
disp('Indices of neurons with significant mutual spatial information excluding track ends:');
disp(significantIndicesSubset);

%%
sigAll    = significantIndices(:);
sigSubset = significantIndicesSubset(:);

% Put into one table. If lengths differ, pad the shorter with NaN so both columns align.
nRows = max(numel(sigAll), numel(sigSubset));
sigAll(end+1:nRows, 1)    = NaN;
sigSubset(end+1:nRows, 1) = NaN;

Tsig = table(sigAll, sigSubset, ...
    'VariableNames', {'significantIndices','significantIndicesSubset'});

% Build output filename from trackingFile
[~, baseName, ~] = fileparts(trackingFile);  % removes .csv (or any extension)
outMatName = baseName + "_cellsWithSignificantSpatialInformation.mat";
outMatPath = fullfile(directoryPath, outMatName);

% Save the table
save(outMatPath, "Tsig");

disp("Saved table to: " + outMatPath);
