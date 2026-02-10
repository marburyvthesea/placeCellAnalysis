%% Paths
baseDir = "/Users/johnmarshall/Documents/Analysis/miniscope_analysis/miniscopeLinearTrack/mIAnalysis/";

csvPath = fullfile(baseDir, "388_run4_20250703_motion_correctedcellTracesAlignedToTracking.csv");
matPath = fullfile(baseDir, "388_run4_20250703_motion_correctedcellTracesAlignedToTracking_cellsWithSignificantSpatialInformation.mat");

%% Load tracking CSV (as table so we can use X_coor)
T = readtable(csvPath);

% If your CSV has an extra header row you previously removed with "fileData(1,:) = [];",
% you can uncomment the next line:
% T(1,:) = [];

% Grab X coordinate
assert(any(strcmp("X_coor", T.Properties.VariableNames)), "Couldn't find column 'X_coor' in the CSV.");
xCoords = T.X_coor;

%% Load significant cell indices table from MAT
S = load(matPath);

% --- Option A: auto-pick the first table found in the MAT file ---
sigTableName = "";
fns = fieldnames(S);
for k = 1:numel(fns)
    if istable(S.(fns{k}))
        sigTableName = fns{k};
        break;
    end
end
assert(sigTableName ~= "", "No table found inside the .mat file. Check what variables are saved in it.");
sigTable = S.(sigTableName);

% Choose which indices column to use: "significantIndices" or "significantIndicesSubset"
whichSigColumn = "significantIndicesSubset";   % or "significantIndices"

assert(ismember(whichSigColumn, sigTable.Properties.VariableNames), ...
    "Column '%s' not found. Available: %s", whichSigColumn, strjoin(sigTable.Properties.VariableNames, ", "));

sigIdx = sigTable.(whichSigColumn);

% Make sure it's a numeric column vector
if iscell(sigIdx), sigIdx = cell2mat(sigIdx); end
sigIdx = sigIdx(:);

% Your sigIdx is MATLAB-style (1 -> cell_0), so NO +1 needed here.
% Only do the +1 conversion if you *know* it's 0-based.
% If you want to keep the heuristic, leave this:
if any(sigIdx == 0)
    sigIdx = sigIdx + 1;
end
%% Identify neuron activity columns in the CSV (robust to missing cell_* labels)
varNames = string(T.Properties.VariableNames);

% Take only columns that start with "cell_"
neuronVarNames = varNames(startsWith(varNames, "cell_"));
assert(~isempty(neuronVarNames), "No neuron columns like 'cell_*' found in the CSV.");

% Sort these columns by the numeric suffix (cell_0, cell_1, cell_2, ...)
cellNums = nan(size(neuronVarNames));
for i = 1:numel(neuronVarNames)
    tok = regexp(neuronVarNames(i), "^cell_(\d+)$", "tokens", "once");
    if ~isempty(tok)
        cellNums(i) = str2double(tok{1});
    end
end

% Keep only properly-formatted names cell_<number>
keep = ~isnan(cellNums);
neuronVarNames = neuronVarNames(keep);
cellNums = cellNums(keep);

[~, ord] = sort(cellNums);
neuronVarNames = neuronVarNames(ord);   % ordered existing neuron columns
% cellNums(ord) tells you which original cell numbers exist

% Now select ONLY the significant neurons by index into this ordered list
assert(all(sigIdx >= 1 & sigIdx <= numel(neuronVarNames)), ...
    "Some sigIdx are out of range. sigIdx must be between 1 and %d (number of existing cell_* columns).", numel(neuronVarNames));

neuronActivity = T{:, neuronVarNames(sigIdx)};


%% Bin along X and compute mean activity per bin per neuron
numBins = 32;

valid = ~isnan(xCoords) & all(~isnan(neuronActivity), 2); % optional: require activity not NaN
x = xCoords(valid);
A = neuronActivity(valid, :);

% Bin edges along observed X range
binEdges = linspace(min(x), max(x), numBins + 1);

binned = NaN(numBins, size(A,2));
for b = 1:numBins
    inBin = x >= binEdges(b) & x < binEdges(b+1);
    if any(inBin)
        binned(b,:) = mean(A(inBin,:), 1, "omitnan");
    end
end

% Optional: normalize each neuron (row-wise in the eventual heatmap) for contrast
% (common for place field plots)
% binned = binned ./ max(binned, [], 1, "omitnan");

%% Sort neurons by peak bin
[~, peakBin] = max(binned, [], 1, "omitnan");
[~, sortOrder] = sort(peakBin);

sortedBinned = binned(:, sortOrder);

% Y labels: use the neuron indices (from your MAT table) in sorted order
sortedNeuronLabels = sigIdx(sortOrder);

%% Plot heatmap
figure;
imagesc(1:numBins, 1:numel(sortedNeuronLabels), sortedBinned');  % neurons on Y, bins on X
set(gca, "YDir", "normal");
xlabel("Spatial bin (along X)");
ylabel("Place cell index");
title("Place cells sorted by peak firing bin");
colorbar;

% Label rows with your cell “names” (indices)
yticks(1:numel(sortedNeuronLabels));
yticklabels(arrayfun(@num2str, sortedNeuronLabels, "UniformOutput", false));

% Colormap (choose one)
% colormap("hot");
colormap("parula");  % MATLAB default
% colormap("jet");

% Optional: tighten contrast
% caxis([min(sortedBinned(:),[],"omitnan"), prctile(sortedBinned(:), 99)]);

%% Make a peak-normalized version (each cell scaled to its own max across bins)
% sortedBinned is numBins x numCells

rowMax = max(sortedBinned, [], 1, "omitnan");   % 1 x numCells (max per cell)
rowMax(rowMax == 0 | isnan(rowMax)) = NaN;      % avoid divide-by-zero

sortedBinned_normToPeak = sortedBinned ./ rowMax;  % implicit expansion (R2016b+)

%% Plot peak-normalized heatmap
figure;

imagesc(1:numBins, 1:numel(sortedNeuronLabels), sortedBinned_normToPeak');  % cells on Y
set(gca, "YDir", "normal");
xlabel("Spatial bin (along X)");
ylabel("Place cell index");
title("Place cells sorted by peak bin (peak-normalized per cell)");
colorbar;

yticks(1:numel(sortedNeuronLabels));
yticklabels(arrayfun(@num2str, sortedNeuronLabels, "UniformOutput", false));

colormap("parula");  % or whatever you like
caxis([0 1]);        % because it's normalized to peak
