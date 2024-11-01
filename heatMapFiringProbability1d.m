
%%load csv file with peaks and tracking data
tableOutput = readtable('/Users/johnmarshall/Documents/MATLAB/PlaceCellAnalysis/233_09_29_sorted_output_98_cells.mat');

%%input list of cell indicies to plot
list_to_plot = [1, 2, 3];

placeCellFiringData = tableOutput{:, list_to_plot};

locationData = tableOutput{:,[101, 99, 100]};

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
colormap('jet'); 
colorbar;

% Add labels and title
xlabel('X Coordinate');
ylabel('Neuron Index');
title('Stacked 1D Heatmaps for All Neurons');

