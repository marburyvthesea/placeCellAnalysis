
%neuron_data =tableOutput{:,[15, 101, 99, 100]};

% Assuming 'neuron_data' is a 11858x4 double with columns: active(0/1), x, y

% Define the bin edges for x and y
x_bin_edges = linspace(min(neuron_data(:, 3)), max(neuron_data(:, 3)), 50); % Adjust the number of bins as needed
y_bin_edges = linspace(min(neuron_data(:, 4)), max(neuron_data(:, 4)), 50); % Adjust the number of bins as needed

% Calculate 2D histogram
histogram_2d = histcounts2(neuron_data(:, 3), neuron_data(:, 4), x_bin_edges, y_bin_edges);

% Calculate the probability of the neuron being active in each bin
probability_matrix = histogram_2d / sum(histogram_2d(:));

% Create a 2D heatmap
figure;
imagesc(x_bin_edges, y_bin_edges, probability_matrix');
colormap('hot');
colorbar;

% Add labels and title
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('2D Heatmap of Neuron Activity Probability');
