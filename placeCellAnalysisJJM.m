
%% converting time row to seconds 

pos(1,1:end)=pos(1,1:end)./100;

%% check if cells have mean firing rate >0.01Hz during periods of movement

%% calculate 2d velocity 
time = pos(1, :); % Time in seconds
x_position = pos(2, :); % X coordinate in cm
y_position = pos(3, :); % Y coordinate in cm
% Calculate velocities using central difference formula
dt = diff(time);
x_velocity = [NaN, diff(x_position) ./ dt, NaN];
y_velocity = [NaN, diff(y_position) ./ dt, NaN];
% Compute overall 2D velocity
velocity_2d = sqrt(x_velocity.^2 + y_velocity.^2);

%% calculate spike rate using 1 second sliding window
% Extract time information
time = pos(1, :);
% Define the 1-second window
window_size = round(1 / mean(diff(time)));
% Create a rectangular window for convolution
window = ones(1, window_size);
% Perform sliding window sum using convolution
event_rate = conv2(spikes', window, 'valid') / window_size;
% Normalize by the width of the window to get events/second
event_rate = event_rate / window_size;
% Display or use 'event_rate' as needed

%% compute mutual information 
%% track position
% bin track into 4cm bins
bin_edges = 0:4:max(x_position)+4;
% Compute histogram counts and bin edges
[counts, edges] = histcounts(x_position, bin_edges);
% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
% Calculate probability of the object occupying each bin
probability_per_bin = counts / sum(counts);


%% firing probability in bins 
% Assuming 'event_rate' is a 98x11839 array and 'x_position' is a 1x11858 array

% Initialize a matrix to store the probabilities for each element
probability_matrix = zeros(size(event_rate, 1), length(counts));
%%
% Loop over each element in 'event_rate'
for i = 1:size(event_rate, 1)
    % get indicies where cell event rate is nonzero 
    nonzero_indices = find(event_rate(i, :) ~= 0);
    
    %get x location at these indicies 
    x_position_when_active=x_position(nonzero_indices);

    % get the counts for when the mouse is active in each bin, in # of
    % samples
    [event_counts_per_x_bin, edges] = histcounts(x_position_when_active, bin_edges);

    % Calculate the probability for each bin
    probability_per_bin = event_counts_per_x_bin ./ length(event_rate(i,:));
    probability_matrix(i, :) = probability_per_bin;
end

% Display or use 'probability_matrix' as needed

%% calculate mutual information per cell


% probability of calcium event over occupany probability
expanded_probability_per_bin = repmat(probability_per_bin, size(probability_matrix, 1), 1);
p_event_over_p_occupancy = probability_matrix ./ expanded_probability_per_bin;

% sum of event probabilities in each bin
MI_percell = zeros(1, size(probability_matrix,1));
for i = 1:size(probability_matrix,1)
    % 
    sum_eventprobability = sum(probability_matrix(1,:));

    sum_occupancyprobability = sum(probability_per_bin); 

    MI_perbin = zeros(1, size(probability_per_bin, 2));
    for j = 1:size(probability_per_bin, 2)

        MI_thisbin = p_event_over_p_occupancy(i,j)*log(p_event_over_p_occupancy(i,j)/(sum_eventprobability+sum_occupancyprobability));
        MI_perbin(1, j) = MI_thisbin; 

    end

    MI_this_cell = sum(MI_perbin) ;

    MI_percell(1, i) = MI_this_cell ; 
end























