

list_to_plot = forward_place_cells_99ci;

placeCellFiringData = tableOutput{:, list_to_plot};

locationData = tableOutput{:,[101, 99, 100]};
%%
% Create a custom colormap from blue to yellow with red in the middle
numColors = 256; % Number of colors in the colormap
halfNumColors = numColors / 2;
% Blue to Red for the first half
blueToRed = [linspace(0, 1, halfNumColors)', zeros(halfNumColors, 1), linspace(1, 0, halfNumColors)'];
% Red to Yellow for the second half
redToYellow = [ones(halfNumColors, 1), linspace(0, 1, halfNumColors)', zeros(halfNumColors, 1)];
% Combine both halves
customColormap = [blueToRed; redToYellow];

% Assuming 'placeCellFiringData' is a 11858x18 double and 'locationData' is a Nx3 double
%% plot 2d heatmaps for each neuron
% Define bin edges for x and y
x_bin_edges = linspace(min(locationData(:, 2)), max(locationData(:, 2)), 50);
y_bin_edges = linspace(min(locationData(:, 3)), max(locationData(:, 3)), 50);

% Iterate through each neuron column in placeCellFiringData
for neuronIndex = 1:size(placeCellFiringData, 2)
    % Extract firing data for the current neuron
    neuronActivity = placeCellFiringData(:, neuronIndex);
    
    % Filter location data based on active neurons
    activeLocations = locationData(logical(neuronActivity), 2:3);
    
    % Calculate 2D histogram for (x, y) coordinates
    histogram_2d = histcounts2(activeLocations(:, 1), activeLocations(:, 2), x_bin_edges, y_bin_edges);
    
    % Calculate the probability of the neuron being active at each (x, y) coordinate
    probability_matrix = histogram_2d / sum(histogram_2d(:));
    
    % Create a heatmap for the current neuron
    figure;
    imagesc(x_bin_edges, y_bin_edges, probability_matrix');
    colormap(customColormap); % Use the custom colormap
    colorbar;
    
    % Add labels and title
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title(['2D Heatmap for Neuron ' num2str(neuronIndex)]);
end

%%
% Define bin edges for x
x_bin_edges = linspace(min(locationData(:, 2)), max(locationData(:, 2)), 40);

% Get the number of neurons
numNeurons = size(placeCellFiringData, 2);

% Initialize a matrix to hold all the heatmaps
allNeuronsHeatmap = zeros(numNeurons, length(x_bin_edges)-1);

% Iterate through each neuron column in placeCellFiringData
for neuronIndex = 1:numNeurons
    % Extract firing data for the current neuron
    neuronActivity = placeCellFiringData(:, neuronIndex);
    
    % Filter x location data based on active neurons
    activeXLocations = locationData(logical(neuronActivity), 2);
    
    % Calculate 1D histogram for x coordinates
    histogram_1d = histcounts(activeXLocations, x_bin_edges);
    
    % Calculate the probability of the neuron being active at each x coordinate
    probability_vector = histogram_1d / sum(histogram_1d);
    
    % Store the probability vector in the matrix
    allNeuronsHeatmap(neuronIndex, :) = probability_vector;
end

% Create a stacked heatmap for all neurons
figure;
imagesc(x_bin_edges(1:end-1), 1:numNeurons, allNeuronsHeatmap);
colormap('jet'); % Or use any other colormap you prefer
colorbar;

% Add labels and title
xlabel('X Coordinate');
ylabel('Neuron Index');
title('Stacked 1D Heatmaps for All Neurons');

