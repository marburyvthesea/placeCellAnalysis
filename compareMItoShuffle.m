%% calc neurons with significant spatial information 
% Define the HDF5 file path and dataset name for shuffled mutual information data
h5ResultsFilePath = strcat(filePath, result{1}, '_MI_results.h5');
datasetName = '/MI_perCellAllShuffles';

% Load the shuffled MI data for all 1000 shuffles
MI_perCellAllShuffles = h5read(h5ResultsFilePath, datasetName);

% Compute the 95% confidence interval for the shuffled data
lowerBound = prctile(MI_perCellAllShuffles, 2.5, 2); % 2.5th percentile
upperBound = prctile(MI_perCellAllShuffles, 97.5, 2); % 97.5th percentile

% Find indices of neurons where actual MI falls outside the 95% CI of shuffled data
significantIndices = find(MI_perCell < lowerBound | MI_perCell > upperBound);

% Display the indices of significant neurons
disp('Indices of neurons with significant mutual spatial information:');
disp(significantIndices);


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
