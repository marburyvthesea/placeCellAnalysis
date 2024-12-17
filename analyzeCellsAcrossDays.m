%%have a data structure for each mouse with place cells/day in
%%significantIndiciesSubset
 
%have "key" that has indicies across days
%use row indicies of cell_registered_struct.cell_to_index_map

%compile "unqiue" global key for each cell that is a place cell across all
%days that stores its identity across days

%[z, numDays] = size(dataStruct);

%%
%sigIndiciesEachDay=zeros(numDays);
%for i = 1:numDays=
%    [2;26;51;93]

%for each cell that was a place cell in any day, compile it's spatial
%information and to firing bins across all days

%%
neuron_matrix = [cell_registered_struct.cell_to_index_map(:,1), ...
    cell_registered_struct.cell_to_index_map(:,2), ...
    cell_registered_struct.cell_to_index_map(:,3)];

%get sig indicies from data structure 
day2_significantIndicesSubset = dataStruct(1).significantIndicesSubset;
day3_significantIndicesSubset = dataStruct(2).significantIndicesSubset;
day4_significantIndicesSubset = dataStruct(3).significantIndicesSubset;

% Get the maximum number of rows among the matrices
maxRows = max([size(day2_significantIndicesSubset, 1), ...
               size(day3_significantIndicesSubset, 1), ...
               size(day4_significantIndicesSubset, 1)]);

% Pad each matrix with zeros to match the maximum number of rows
padded_day2 = padarray(day2_significantIndicesSubset, ...
                       [maxRows - size(day2_significantIndicesSubset, 1), 0], ...
                       0, 'post');

padded_day3 = padarray(day3_significantIndicesSubset, ...
                       [maxRows - size(day3_significantIndicesSubset, 1), 0], ...
                       0, 'post');

padded_day4 = padarray(day4_significantIndicesSubset, ...
                       [maxRows - size(day4_significantIndicesSubset, 1), 0], ...
                       0, 'post');
%%
% Horizontally concatenate the padded matrices
concatenatedSigNeurons = [padded_day2, padded_day3, padded_day4];

[numSignificantNeurons, numDays] = size(concatenatedSigNeurons); 

sigNeuronsRegisteredAcrossDays = [];
for day=1:numDays

    for row=1:numSignificantNeurons
        
        thisNeuronThisDay = concatenatedSigNeurons(row, day); 
        
        if thisNeuronThisDay > 0
            rowIndex = find(neuron_matrix(:, day) == thisNeuronThisDay, day);
            thisNeuronAllDays=neuron_matrix(rowIndex,:); 
            sigNeuronsRegisteredAcrossDays = [sigNeuronsRegisteredAcrossDays; thisNeuronAllDays];
        end
    end
end

% Find rows where there aren't any zeros
rowsSigNeuronsInAllDays = all(sigNeuronsRegisteredAcrossDays ~= 0, 2);
% Extract the rows without zeros
sigNeuronsInAllDays = sigNeuronsRegisteredAcrossDays(rowsSigNeuronsInAllDays, :);

%% plot heatmaps 
%% extract calcium info and mouse position for "place cells" aligned across multiple days
pathToDirectoryWithTrackingFiles='/Users/johnmarshall/Documents/Analysis/miniscope_lineartrack/m3_AlignedToTracking_toalign/';
% cell firing information 
day2_dataaligned = dataStruct(3).fileName;
day3_dataaligned = dataStruct(2).fileName;
day4_dataaligned = dataStruct(1).fileName;

concatenatedFiles = [day2_dataaligned; day3_dataaligned; ...
    day4_dataaligned];
[numFiles, ~] = size(concatenatedFiles); 

extractedCalciumSignals={};
extractedPositionInformation={};

for i = 1:numFiles
    disp(i);
    % Load the current file
    fileName = concatenatedFiles(i,:);
    disp(fileName);
    significantIndicesSubset = sigNeuronsInAllDays(:,i) + 1;  % Add 1 to adjust for indexing
    significantIndicesOrigVariableNames = sigNeuronsInAllDays(:,i);

    % Read the .csv file into a table or array
    fileDataTable = readtable(strcat(pathToDirectoryWithTrackingFiles, fileName), 'VariableNamesLine', 1);
    coordinateColumns = {'X_coor', 'Y_coor'};
    coordinateTable = fileDataTable(2:end, coordinateColumns);
    fileData = readmatrix(strcat(pathToDirectoryWithTrackingFiles, fileName));  % Using readmatrix to handle .csv format directly
    %remove 1st row (variable names) and 1st column (timestamps)
    fileData(1, :) = [];
    fileData(:, 1) = [];

    % Extract the specified columns (significant cells)
    significantData = fileData(:, significantIndicesSubset);
    significantIndicesOrigVariableNames = cellstr(string(significantIndicesOrigVariableNames));
    significantIndicesOrigVariableNames = matlab.lang.makeUniqueStrings(significantIndicesOrigVariableNames);
    significantTable = array2table(significantData, 'VariableNames', significantIndicesOrigVariableNames);

    % Store the extracted data in the cell array
    extractedCalciumSignals{i} = significantTable;
    extractedPositionInformation{i} = coordinateTable;
end

%% bin neural activity on different days

numBins = 32;
% Initialize a cell array to store the averaged neural activity for each file
binnedActivity = cell(numFiles, 1);
for i = 1:numFiles
    xCoords = extractedPositionInformation{i}.X_coor;
    neuronActivity = extractedCalciumSignals{i}; 

    % Define bin edges based on the maximum X coordinate in the current file
    maxX = max(xCoords);
    binEdges = linspace(0, maxX, numBins + 1);

    % Create an empty table to store binned activity with matching column names
    binnedTable = array2table(NaN(numBins, width(neuronActivity)), ...
                              'VariableNames', neuronActivity.Properties.VariableNames);

    % Loop over each spatial bin
    for bin = 1:numBins
        % Find the indices of xCoords that fall within the current bin
        inBinIdx = xCoords >= binEdges(bin) & xCoords < binEdges(bin + 1);

        % Calculate the average activity for each neuron within this bin
        if any(inBinIdx)
            % Use brace indexing to extract numeric data and compute mean
            binnedTable{bin, :} = mean(neuronActivity{inBinIdx, :}, 1);
        end
    end

    % Store the binned table in the cell array
    binnedActivity{i} = binnedTable;

end

%%
% Number of days (or entries in binnedActivity)
numDays = length(binnedActivity);

% Extract the number of columns from the first table in binnedActivity
numColumns = width(binnedActivity{1});

% Iterate over each column
for colIdx = 1:numColumns
    % Create a figure for the heatmaps for the current column
    figure;
    tiledlayout(1, numDays, 'Padding', 'compact', 'TileSpacing', 'compact'); % Arrange side-by-side

    for dayIdx = 1:numDays
        % Extract the data for the current column and day
        dataForHeatmap = binnedActivity{dayIdx}{:, colIdx}; % Extract column as numeric vector

        % Create a subplot for the current day
        nexttile;
        imagesc(dataForHeatmap); % Generate heatmap
        colormap('hot'); % Use 'hot' colormap
        colorbar; % Add a colorbar
        caxis([min(dataForHeatmap(:)), max(dataForHeatmap(:))]); % Set consistent color limits

        % Add labels and title
        title(['Day ' num2str(dayIdx)]);
        xlabel('Spatial Bin');
        ylabel('Magnitude');
    end

    % Add a title for the column
    sgtitle(['Heatmaps for Column ' binnedActivity{1}.Properties.VariableNames{colIdx}]);
end









