%%heat map from MI information dataStructure
%%add neurons registered across all days 

%%
pathToDirectoryWithTrackingFiles = '/Users/johnmarshall/Documents/Analysis/miniscope_lineartrack/m3_AlignedToTracking_toalign/' ; 

% Initialize a cell array to store the extracted data with padding for each file
extractedData = cell(length(dataStruct), 1);
extractedDataTables = cell(length(dataStruct), 1); 

% Loop over each file in dataStruct
for i = 1:length(dataStruct)
    % Load the current file
    fileName = dataStruct(i).fileName;
    significantIndicesSubset = dataStruct(i).significantIndicesSubset + 1;  % Add 1 to adjust for indexing

    % Read the .csv file into a table or array
    fileDataTable = readtable(strcat(pathToDirectoryWithTrackingFiles, fileName), 'VariableNamesLine', 1);
    fileData = readmatrix(strcat(pathToDirectoryWithTrackingFiles, fileName));  % Using readmatrix to handle .csv format directly
    fileData(1, :) = [];

    % Extract the specified columns (significant cells)
    significantData = fileData(:, significantIndicesSubset);

    % Store the extracted data in the cell array
    extractedData{i} = significantData;
    extractedDataTables{i}  = fileDataTable;
end

%% Determine the maximum number of rows and columns across all files
maxRows = max(cellfun(@(x) size(x, 1), extractedData));
maxCols = max(cellfun(@(x) size(x, 2), extractedData));

% Pad each extracted dataset to match the maximum number of rows and columns
for i = 1:length(extractedData)
    data = extractedData{i};
    % Pad rows if needed
    if size(data, 1) < maxRows
        data = [data; NaN(maxRows - size(data, 1), size(data, 2))];
    end
    % Pad columns if needed
    if size(data, 2) < maxCols
        data = [data, NaN(size(data, 1), maxCols - size(data, 2))];
    end
    % Store the padded data back in the cell array
    extractedData{i} = data;
end

% Combine the data into a single 3D array (maxRows x maxCols x numFiles)
numFiles = length(dataStruct);
paddedDataArray = NaN(maxRows, maxCols, numFiles);  % Initialize a 3D array with NaNs

for i = 1:numFiles
    paddedDataArray(:, :, i) = extractedData{i};  % Store each file's data in the 3D array
end

%% extract location information for each file 

% Initialize a cell array to store the extracted coordinates with padding for each file
extractedCoordinates = cell(length(dataStruct), 1);

% Loop over each file in dataStruct to extract the .X_coor and .Y_coor columns
for i = 1:length(dataStruct)
    % Load the current file
    fileName = dataStruct(i).fileName;

    % Read the .csv file into a table to access column names
    %fileData = readtable(strcat(pathToDirectoryWithTrackingFiles, fileName)); 
    fileData = readtable(strcat(pathToDirectoryWithTrackingFiles, fileName), 'VariableNamesLine', 1);
    fileData(1, :) = [];
    
    % Extract .X_coor and .Y_coor columns if they exist
    coordinatesData = fileData{:, {'X_coor', 'Y_coor'}};

    % Store the extracted coordinates in the cell array
    extractedCoordinates{i} = coordinatesData;
end

% Determine the maximum number of rows across all files (for padding)
maxRows = max(cellfun(@(x) size(x, 1), extractedCoordinates));
numFiles = length(dataStruct);

% Pad each extracted dataset to match the maximum number of rows
for i = 1:numFiles
    data = extractedCoordinates{i};
    if size(data, 1) < maxRows
        % Pad with NaNs to match the maximum row length
        data = [data; NaN(maxRows - size(data, 1), size(data, 2))];
    end
    % Store the padded data back in the cell array
    extractedCoordinates{i} = data;
end

% Combine the data into a single 3D array (maxRows x 2 x numFiles)
paddedCoordinatesArray = NaN(maxRows, 2, numFiles);  % Initialize a 3D array with NaNs

for i = 1:numFiles
    paddedCoordinatesArray(:, :, i) = extractedCoordinates{i};  % Store each file's data in the 3D array
end

%% plot heatmaps 
% Define the number of spatial bins
numBins = 32;

% Initialize a cell array to store the averaged neural activity for each file
binnedActivity = cell(length(dataStruct), 1);

% Loop through each file's data in paddedDataArray and paddedCoordinatesArray
for i = 1:length(dataStruct)
    % Extract the X coordinates and neuron activity for the current file
    xCoords = paddedCoordinatesArray(:, 1, i);
    neuronActivity = paddedDataArray(:, :, i);

    % Remove rows with NaNs in xCoords (padding rows)
    validIdx = ~isnan(xCoords);
    xCoords = xCoords(validIdx);
    neuronActivity = neuronActivity(validIdx, :);

    % Define bin edges based on the maximum X coordinate in the current file
    maxX = max(xCoords);
    binEdges = linspace(0, maxX, numBins + 1);

    % Initialize an array to store the binned average activity for each neuron
    binnedActivity{i} = NaN(numBins, size(neuronActivity, 2));

    % Loop over each spatial bin
    for bin = 1:numBins
        % Find the indices of xCoords that fall within the current bin
        inBinIdx = xCoords >= binEdges(bin) & xCoords < binEdges(bin + 1);

        % Calculate the average activity for each neuron within this bin
        if any(inBinIdx)
            binnedActivity{i}(bin, :) = mean(neuronActivity(inBinIdx, :), 1);
        end
    end
end


%% Combine all neurons into a single matrix for plotting
combinedBinnedActivity = [];
combinedNeuronIndices = [];

for i = 1:length(binnedActivity)
    % Extract the binned activity for the current file
    currentDayActivity = binnedActivity{i};
    currentNeuronIndices = dataStruct(i).significantIndicesSubset;  % Neuron indices for the current file

    % Find valid columns (neurons) that are not entirely NaN
    validColumns = any(~isnan(currentDayActivity), 1);
    currentDayActivity = currentDayActivity(:, validColumns);
    currentNeuronIndices = currentNeuronIndices(validColumns);  % Keep only valid indices

    % Convert neuron indices to strings and make them unique
    neuronLabels = matlab.lang.makeUniqueStrings(string(currentNeuronIndices), combinedNeuronIndices);

    % Append current day's data and unique indices to the combined matrices
    combinedBinnedActivity = [combinedBinnedActivity, currentDayActivity];
    combinedNeuronIndices = [combinedNeuronIndices; neuronLabels(:)];
end

% Find the index of the maximum firing bin for each neuron
[~, maxBinIndices] = max(combinedBinnedActivity, [], 1);

% Sort neurons based on the location of their maximum firing bin
[~, sortOrder] = sort(maxBinIndices);

% Reorder combinedBinnedActivity and neuron indices by the sorted order of neurons
sortedBinnedActivity = combinedBinnedActivity(:, sortOrder);
sortedNeuronIndices = combinedNeuronIndices(sortOrder);

%% Plot the heatmap
figure;
imagesc(sortedBinnedActivity');  % Transpose to have neurons on the y-axis and bins on the x-axis
colorbar;

% Define a custom colormap going from blue to green to yellow to red
numColors = 256;  % Number of color steps
customColormap = [linspace(0, 1, numColors)', linspace(0, 1, numColors)', linspace(1, 0, numColors)'];

% Apply the custom colormap
colormap(customColormap);

% Adjust color limits to improve contrast (optional)
caxis([min(sortedBinnedActivity(:)), max(sortedBinnedActivity(:)) * 0.7]);  % Adjust the contrast

% Add y-axis labels with unique neuron indices
yticks(1:size(sortedBinnedActivity, 2));  % Set y-ticks for all rows
yticklabels(cellstr(sortedNeuronIndices));  % Convert unique neuron indices to cell array of strings

% Add plot labels and title
title('Heatmap of Binned Neural Activity with Unique Neuron Indices');
xlabel('Spatial Bin');
ylabel('Neuron Index');
set(gca, 'YDir', 'normal');  % Correct the y-axis direction

% Save the heatmap as a .tif file
%saveas(gcf, 'heatmap_binned_neural_activity_.svg');
%fileName = 'heatmap_binned_neural_activity_.tif';
%exportgraphics(gcf, fileName, 'Resolution', 300); 

%%plot showing boundaries by day
%% Combine all neurons into a single matrix for plotting
combinedBinnedActivity = [];
combinedNeuronIndices = [];  % To store neuron indices for the y-axis labels

% Track the start and end indices for each day's data
dayBoundaries = zeros(1, numel(binnedActivity) + 1);
currentIndex = 1;

sortedDayActivities = cell(1, numel(binnedActivity)); % Store sorted activities for each day

for i = 1:numel(binnedActivity)
    % Extract the binned activity for the current day
    currentDayActivity = binnedActivity{i};
    currentNeuronIndices = dataStruct(i).significantIndicesSubset;  % Neuron indices for the current day

    % Find valid columns (neurons) that are not entirely NaN
    validColumns = any(~isnan(currentDayActivity), 1);
    currentDayActivity = currentDayActivity(:, validColumns);
    currentNeuronIndices = currentNeuronIndices(validColumns);  % Keep only valid indices

    % Arrange neurons by the location of their maximum firing bin
    [~, maxBinIndices] = max(currentDayActivity(:,:), [], 1, 'omitnan');
    [~, sortOrder] = sort(maxBinIndices);

    % Sort neurons for the current day
    sortedDayActivity = currentDayActivity(:, sortOrder);
    sortedNeuronIndices = currentNeuronIndices(sortOrder);

    % Store the sorted activity for this day
    sortedDayActivities{i} = sortedDayActivity;

    % Append current day's data to the combined matrix
    combinedBinnedActivity = [combinedBinnedActivity, sortedDayActivity(:,:)];
    combinedNeuronIndices = [combinedNeuronIndices; sortedNeuronIndices];

    % Record the boundary for this day
    currentIndex = currentIndex + size(sortedDayActivity, 2);
    dayBoundaries(i + 1) = currentIndex - 1;
end

%% Plot the heatmap
figure;
imagesc(combinedBinnedActivity');  % Transpose to have neurons on the y-axis and bins on the x-axis
colorbar;

% Define a custom colormap going from blue to green to yellow to red
numColors = 256;  % Number of color steps
customColormap = [linspace(0, 1, numColors)', linspace(0, 1, numColors)', linspace(1, 0, numColors)'];

% Apply the custom colormap
colormap(customColormap);

% Adjust color limits to improve contrast (optional)
caxis([min(combinedBinnedActivity(:)), max(combinedBinnedActivity(:)) * 0.7]);  % Adjust the contrast

% Add black lines to separate data from different days
hold on;
for i = 2:numel(dayBoundaries) - 1
    yLine = dayBoundaries(i) + 0.5;  % Add 0.5 for the center of the gap
    line([0.5, size(combinedBinnedActivity, 1) + 0.5], [yLine, yLine], 'Color', 'k', 'LineWidth', 2);
end
hold off;

% Add y-axis labels with neuron indices
yticks(1:size(combinedBinnedActivity, 2));  % Set y-ticks for all rows
yticklabels(arrayfun(@num2str, combinedNeuronIndices, 'UniformOutput', false));  % Convert neuron indices to strings

% Add plot labels and title
title('Heatmap of Binned Neural Activity with Neuron Indices');
xlabel('Spatial Bin');
ylabel('Neuron Index');
set(gca, 'YDir', 'normal');  % Correct the y-axis direction


%%
% Loop through the columns of sigNeuronsInAllDays and assign to the structure
for colIdx = 1:size(sigNeuronsInAllDays, 2)
    dataStruct(colIdx).sigNeuronsInAllDays = sigNeuronsInAllDays(:, colIdx);
end
%%
% Initialize a cell array to store the extracted data with padding for each file
extractedDataAllDays = cell(length(dataStruct), 1);
extractedDataTablesAllDays = cell(length(dataStruct), 1); 

% Loop over each file in dataStruct
for i = 1:length(dataStruct)
    % Load the current file
    fileName = dataStruct(i).fileName;
    significantIndicesAllDaysSubset = dataStruct(i).sigNeuronsInAllDays + 1;  % Add 1 to adjust for indexing

    % Read the .csv file into a table or array
    fileDataTable = readtable(strcat(pathToDirectoryWithTrackingFiles, fileName), 'VariableNamesLine', 1);
    fileData = readmatrix(strcat(pathToDirectoryWithTrackingFiles, fileName));  % Using readmatrix to handle .csv format directly
    fileData(1, :) = [];

    % Extract the specified columns (significant cells)
    significantData = fileData(:, significantIndicesAllDaysSubset);

    % Store the extracted data in the cell array
    extractedDataAllDays{i} = significantData;
    extractedDataTablesAllDays{i} = fileDataTable; 
end

%% Determine the maximum number of rows and columns across all files
maxRowsAllDays = max(cellfun(@(x) size(x, 1), extractedDataAllDays));
maxColsAllDays = max(cellfun(@(x) size(x, 2), extractedDataAllDays));

% Pad each extracted dataset to match the maximum number of rows and columns
for i = 1:length(extractedDataAllDays)
    data = extractedDataAllDays{i};
    % Pad rows if needed
    if size(data, 1) < maxRowsAllDays
        data = [data; NaN(maxRowsAllDays - size(data, 1), size(data, 2))];
    end
    % Pad columns if needed
    if size(data, 2) < maxColsAllDays
        data = [data, NaN(size(data, 1), maxColsAllDays - size(data, 2))];
    end
    % Store the padded data back in the cell array
    extractedDataAllDays{i} = data;
end

% Combine the data into a single 3D array (maxRows x maxCols x numFiles)
numFiles = length(dataStruct);
paddedDataArrayAllDays = NaN(maxRowsAllDays, maxColsAllDays, numFiles);  % Initialize a 3D array with NaNs

for i = 1:numFiles
    paddedDataArrayAllDays(:, :, i) = extractedDataAllDays{i};  % Store each file's data in the 3D array
end

%% extract location information for each file
% should be able to load previous coordinate information from above 

%% plot heatmaps 
% Define the number of spatial bins
numBins = 32;

% Initialize a cell array to store the averaged neural activity for each file
binnedActivityAllDays = cell(length(dataStruct), 1);

% Loop through each file's data in paddedDataArray and paddedCoordinatesArray
for i = 1:length(dataStruct)
    % Extract the X coordinates and neuron activity for the current file
    xCoords = paddedCoordinatesArray(:, 1, i);
    neuronActivity = paddedDataArrayAllDays(:, :, i);

    % Remove rows with NaNs in xCoords (padding rows)
    validIdx = ~isnan(xCoords);
    xCoords = xCoords(validIdx);
    neuronActivity = neuronActivity(validIdx, :);

    % Define bin edges based on the maximum X coordinate in the current file
    maxX = max(xCoords);
    binEdges = linspace(0, maxX, numBins + 1);

    % Initialize an array to store the binned average activity for each neuron
    binnedActivityAllDays{i} = NaN(numBins, size(neuronActivity, 2));

    % Loop over each spatial bin
    for bin = 1:numBins
        % Find the indices of xCoords that fall within the current bin
        inBinIdx = xCoords >= binEdges(bin) & xCoords < binEdges(bin + 1);

        % Calculate the average activity for each neuron within this bin
        if any(inBinIdx)
            binnedActivityAllDays{i}(bin, :) = mean(neuronActivity(inBinIdx, :), 1);
        end
    end
end

%%
%% Combine all neurons into a single matrix for plotting across all days
% Concatenate all columns (neurons) from each file's binned activity
combinedBinnedActivityAllDays = cat(2, binnedActivityAllDays{:});
% Find columns that are not entirely NaN
validColumns = any(~isnan(combinedBinnedActivityAllDays), 1);
% Keep only the columns that contain at least some data
combinedBinnedActivityAllDays = combinedBinnedActivityAllDays(:, validColumns);

%Arrange cells by the location of their maximum firing bin
% Find the index of the maximum firing bin for each neuron
[~, maxBinIndices] = max(combinedBinnedActivityAllDays, [], 1);

% Sort neurons based on the location of their maximum firing bin
[~, sortOrder] = sort(maxBinIndices);

% Reorder combinedBinnedActivity by the sorted order of neurons
sortedBinnedActivityAllDays = combinedBinnedActivityAllDays(:, sortOrder);

%% Plot the heatmap
figure;
imagesc(sortedBinnedActivityAllDays');  % Transpose to have neurons on the y-axis and bins on the x-axis
colorbar;

% Define a custom colormap going from blue to green to yellow to red
numColors = 256;  % Number of color steps
customColormap = [linspace(0, 1, numColors)', linspace(0, 1, numColors)', linspace(1, 0, numColors)'];

% Apply the custom colormap
colormap(customColormap);
%colormap('jet');

% Adjust color limits to improve contrast (optional)
caxis([min(sortedBinnedActivityAllDays(:)), max(sortedBinnedActivityAllDays(:)) * 0.7]);  % Adjust the contrast

% Add plot labels and title
title('Heatmap of Binned Neural Activity (Blue-Green-Yellow-Red)');
xlabel('Spatial Bin');
ylabel('Neuron (sorted by maximum firing location)');
set(gca, 'YDir', 'normal');  % Correct the y-axis direction

% Save the heatmap as a .tif file
%saveas(gcf, 'heatmap_binned_neural_activity_.svg');
%fileName = 'heatmap_binned_neural_activity_.tif';
%exportgraphics(gcf, fileName, 'Resolution', 300); 

%%plot showing boundaries by day


%% Plot Place Cells Found in All Days by max firing bin: Combine all neurons into a single matrix for plotting
combinedBinnedActivityAllDays = [];
combinedNeuronIndicesAllDays = [];  % To store the indices of the neurons for the y-axis labels
globalNeuronsAcrossDays = []; 

% Track the start and end indices for each day's data
dayBoundaries = zeros(1, numel(binnedActivityAllDays) + 1);
currentIndex = 1;

sortedDayActivitiesAllDays = cell(1, numel(binnedActivityAllDays)); % Store sorted activities for each day

for i = 1:numel(binnedActivityAllDays)
    % Extract the binned activity for the current day
    currentDayActivity = binnedActivityAllDays{i};
    currentNeuronIndices = dataStruct(i).sigNeuronsInAllDays;  % Indices for the current day

    % Find valid columns (neurons) that are not entirely NaN
    validColumns = any(~isnan(currentDayActivity(:,:)), 1);
    currentDayActivity = currentDayActivity(:, validColumns);
    currentNeuronIndices = currentNeuronIndices(validColumns);  % Keep only valid indices

    % Track neuron positions within the current dataStruct entry
    neuronPositions = find(validColumns);  % Position within original list

    % Arrange neurons by the location of their maximum firing bin
    [~, maxBinIndices] = max(currentDayActivity(:,:), [], 1, 'omitnan');
    [~, sortOrder] = sort(maxBinIndices);

    % Sort neurons for the current day
    sortedDayActivity = currentDayActivity(:, sortOrder);
    sortedNeuronIndices = currentNeuronIndices(sortOrder);
    sortedNeuronPositions = neuronPositions(sortOrder);

    % Store the sorted activity for this day
    sortedDayActivitiesAllDays{i} = sortedDayActivity;

    % Append current day's data to the combined matrix
    combinedBinnedActivityAllDays = [combinedBinnedActivityAllDays, sortedDayActivity(:,:)];
    combinedNeuronIndicesAllDays = [combinedNeuronIndicesAllDays; sortedNeuronIndices(:)];
    globalNeuronsAcrossDays = [globalNeuronsAcrossDays; sortedNeuronPositions(:)];

    % Record the boundary for this day
    currentIndex = currentIndex + size(sortedDayActivity, 2);
    dayBoundaries(i + 1) = currentIndex - 1;
end


%% Plot the heatmap
figure;
fig = gcf;
fig.Position(3) = fig.Position(3) + 200;
ax = gca;  % Get the current axes handle
ax.Position(3) = ax.Position(3) * 0.9;
imagesc(combinedBinnedActivityAllDays');  % Transpose to have neurons on the y-axis and bins on the x-axis
cb = colorbar;
cb.Position(1) = cb.Position(1) + 0.1;

% Define a custom colormap going from blue to green to yellow to red
numColors = 256;  % Number of color steps
customColormap = [linspace(0, 1, numColors)', linspace(0, 1, numColors)', linspace(1, 0, numColors)'];

% Apply the custom colormap
colormap(customColormap);

% Adjust color limits to improve contrast (optional)
caxis([min(combinedBinnedActivityAllDays(:)), max(combinedBinnedActivityAllDays(:)) * 0.7]);  % Adjust the contrast

% Add black lines to separate data from different days
hold on;
for i = 2:numel(dayBoundaries) - 1
    yLine = dayBoundaries(i) + 0.5;  % Add 0.5 for the center of the gap
    line([0.5, size(combinedBinnedActivityAllDays, 1) + 0.5], [yLine, yLine], 'Color', 'k', 'LineWidth', 2);
end
hold off;

% Add y-axis labels with neuron indices
yticks(1:size(combinedBinnedActivityAllDays, 2));  % Set y-ticks for all rows
yticklabels(arrayfun(@num2str, combinedNeuronIndicesAllDays, 'UniformOutput', false));  % Convert neuron indices to strings

% Add secondary y-axis labels (right) with neuron positions
% Annotate each row with its position
for y = 1:length(globalNeuronsAcrossDays)
    text(size(combinedBinnedActivityAllDays, 1) + 1, y, num2str(globalNeuronsAcrossDays(y)), ...
         'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'left');
end


% Add plot labels and title
title('Heatmap of Binned Neural Activity with Neuron Indices');
xlabel('Spatial Bin');
ylabel('Neuron Index (Within Days)');
set(gca, 'YDir', 'normal');  % Correct the y-axis direction

%% Plot Place Cells Found in All Days by max firing bin on day 1, preserve order:
%% Combine all neurons into a single matrix for plotting
combinedBinnedActivityAllDays = [];
combinedNeuronIndicesAllDays = [];  % To store the indices of the neurons for the y-axis labels
globalNeuronsAcrossDays = []; 

% Track the start and end indices for each day's data
dayBoundaries = zeros(1, numel(binnedActivityAllDays) + 1);
currentIndex = 1;

sortedDayActivitiesAllDays = cell(1, numel(binnedActivityAllDays)); % Store sorted activities for each day

% Step 1: Sort neurons in Day 1 by their maximum firing bin
day1Activity = binnedActivityAllDays{1};
validColumnsDay1 = any(~isnan(day1Activity), 1);  % Find valid columns
day1Activity = day1Activity(:, validColumnsDay1);
day1NeuronIndices = dataStruct(1).sigNeuronsInAllDays(validColumnsDay1);

% Determine sorting order based on maximum firing bin in Day 1
[~, maxBinIndicesDay1] = max(day1Activity, [], 1, 'omitnan');
[~, sortOrderDay1] = sort(maxBinIndicesDay1);

% Sort Day 1 activity and neuron indices
sortedDay1Activity = day1Activity(:, sortOrderDay1);
sortedNeuronIndicesDay1 = day1NeuronIndices(sortOrderDay1);

% Append Day 1 data to combined matrices
combinedBinnedActivityAllDays = [combinedBinnedActivityAllDays, sortedDay1Activity];
combinedNeuronIndicesAllDays = [combinedNeuronIndicesAllDays; sortedNeuronIndicesDay1];
globalNeuronsAcrossDays = [globalNeuronsAcrossDays; (1:length(sortOrderDay1))'];  % Track positions

% Step 2: Use Day 1 sorting order to arrange subsequent days
for i = 2:numel(binnedActivityAllDays)
    % Extract activity for the current day
    currentDayActivity = binnedActivityAllDays{i};
    currentNeuronIndices = dataStruct(i).sigNeuronsInAllDays;

    % Find valid columns and keep only neurons matching Day 1 order
    %validColumns = any(~isnan(currentDayActivity), 1);
    %currentDayActivity = currentDayActivity(:, validColumns);
    %currentNeuronIndices = currentNeuronIndices(validColumns);

    % Match the neurons to Day 1 order
    [~, idxInDay1] = ismember(sortedNeuronIndicesDay1, currentNeuronIndices);
    sortedDayActivity = NaN(size(currentDayActivity, 1), length(sortOrderDay1));
    for j = 1:length(idxInDay1)
        if idxInDay1(j) > 0
            sortedDayActivity(:, j) = currentDayActivity(:, idxInDay1(j));
        end
    end

    % Append sorted activity for the current day
    combinedBinnedActivityAllDays = [combinedBinnedActivityAllDays, sortedDayActivity];
    combinedNeuronIndicesAllDays = [combinedNeuronIndicesAllDays; sortedNeuronIndicesDay1];
    globalNeuronsAcrossDays = [globalNeuronsAcrossDays; (1:length(sortOrderDay1))'];  % Reuse positions

    % Record the boundary for this day
    currentIndex = currentIndex + size(sortedDayActivity, 2);
    dayBoundaries(i + 1) = currentIndex - 1;
end

%% Plot the heatmap
figure;
fig = gcf;
fig.Position(3) = fig.Position(3) + 200;  % Expand figure width
ax = gca;
ax.Position(3) = ax.Position(3) * 0.9;   % Adjust axes width

imagesc(combinedBinnedActivityAllDays');  % Transpose for y-axis neurons
cb = colorbar;
cb.Position(1) = cb.Position(1) + 0.05;  % Shift colorbar to the right

% Define a custom colormap
numColors = 256;
customColormap = [linspace(0, 1, numColors)', linspace(0, 1, numColors)', linspace(1, 0, numColors)'];
colormap(customColormap);

% Adjust color limits for contrast
caxis([min(combinedBinnedActivityAllDays(:)), max(combinedBinnedActivityAllDays(:)) * 0.7]);

% Add black lines to separate data from different days
hold on;
for i = 2:numel(dayBoundaries) - 1
    yLine = dayBoundaries(i) + 0.5;
    line([0.5, size(combinedBinnedActivityAllDays, 1) + 0.5], [yLine, yLine], 'Color', 'k', 'LineWidth', 2);
end
hold off;

% Add y-axis labels (neuron indices)
yticks(1:size(combinedBinnedActivityAllDays, 2));
yticklabels(arrayfun(@num2str, combinedNeuronIndicesAllDays, 'UniformOutput', false));

% Add secondary y-axis labels (right side with Day 1 positions)
ax2 = axes('Position', ax.Position, 'YAxisLocation', 'right', 'Color', 'none');
ax2.YLim = ax.YLim;
yticks(1:length(globalNeuronsAcrossDays));
yticklabels(arrayfun(@num2str, globalNeuronsAcrossDays, 'UniformOutput', false));
ax2.XColor = 'none';
ax2.YColor = 'k';

% Add plot labels and title
title('Heatmap of Binned Neural Activity Sorted by Day 1');
xlabel('Spatial Bin');
ylabel('Neuron Index (Within Days)');
set(gca, 'YDir', 'normal');  % Correct y-axis direction




