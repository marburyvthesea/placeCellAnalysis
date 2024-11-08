%%load csv file with peaks and tracking data
run("loadSortedSession.m")

%%
CNMFE_alignedToTracking = dataTable;

%%
frameTimes= dataTable(:,1);
%% converting time row to seconds
timeStrings = frameTimes; % Extract the time strings from the table
timeDurations = duration(extractAfter(timeStrings(:,1).Var1, ' days '), 'Format', 'hh:mm:ss.SSSSSS');
pos = seconds(timeDurations);

cellTraces= dataTable(:,2:189);
X_coor= dataTable.X_coor;
Y_coor= dataTable.Y_coor;


%% check if cells have mean firing rate >0.01Hz during periods of movement

%% calculate 2d velocity 
time = pos(:, 1); % Time in seconds
x_position = X_coor; % X coordinate in cm
y_position = Y_coor; % Y coordinate in cm
% Calculate velocities using central difference formula
dt = diff(time);
x_velocity = [NaN; diff(x_position) ./ dt];
y_velocity = [NaN; diff(y_position) ./ dt];
% Compute overall 2D velocity
velocity_2d = sqrt(x_velocity.^2 + y_velocity.^2);
velocity_2d = [velocity_2d; NaN];

%% calculate spike rate using 1 second sliding window
% label peaks exceeding 2.5 SD threshold
[signalPeaks] = computeSignalPeaks(table2array(cellTraces), 'doMovAvg', 0, 'reportMidpoint', 1, 'numStdsForThresh', 2.5);
%% Extract time information
spikes=signalPeaks;
time = pos;
% Define the 1-second window
window_size = round(1 / mean(diff(time)));
% Create a 1D convolution window to sum across the time dimension
window = ones(window_size, 1);
% Perform sliding window sum along the time dimension using convolution
event_rate = conv2(spikes, window, 'same') / window_size;

%% compute mutual information 
%% track position
% bin track into 4cm bins
bin_edges = 0:4:max(x_position)+4;
% Compute histogram counts and bin edges
[counts, edges] = histcounts(x_position, bin_edges);
% Calculate bin centers
bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
% Calculate probability of the object occupying each bin
probabilityOfMouseOccupyingBin = counts / sum(counts);


%% firing probability in bins 
% Assuming 'event_rate' is a 98x11839 array and 'x_position' is a 1x11858 array
event_rate_t = event_rate'; 
x_position_t = x_position';

% Initialize a matrix to store the probabilities for each element
cellFiringProbabilityPerBin = zeros(size(event_rate_t, 1), length(counts));
% Loop over each element in 'event_rate'
for i = 1:size(event_rate_t, 1)
    % get indicies where cell event rate is nonzero 
    nonzero_indices = find(event_rate_t(i, :) ~= 0);
    
    %get x location at these indicies 
    x_position_when_active=x_position_t(nonzero_indices);

    % get the counts for when the mouse is active in each bin, in # of
    % samples
    [event_counts_per_x_bin, edges] = histcounts(x_position_when_active, bin_edges);

    % Calculate the probability for each bin
    probability_per_bin = event_counts_per_x_bin ./ length(event_rate_t(i,:));
    cellFiringProbabilityPerBin(i, :) = probability_per_bin;
end

% Display or use 'probability_matrix' as needed

%% calculate mutual information per cell

% overall firing probability for each neuron 
neuronFiringProbability = mean(signalPeaks, 1);

% probability of calcium event over occupany probability
expanded_probability_per_bin = repmat(probabilityOfMouseOccupyingBin, size(cellFiringProbabilityPerBin, 1), 1);
p_event_over_p_occupancy = cellFiringProbabilityPerBin ./ expanded_probability_per_bin;

% sum of event probabilities in each bin
MI_perCellperBin = zeros(size(cellFiringProbabilityPerBin,1), size(cellFiringProbabilityPerBin,2));
MI_perCell = zeros(1, size(cellFiringProbabilityPerBin,1));

% for each cell
size(cellFiringProbabilityPerBin,1)
for i = 1:size(cellFiringProbabilityPerBin,1)

    MI_perbin = zeros(1, size(probabilityOfMouseOccupyingBin, 2));
    
    %for each bin  
    for j = 1:size(probabilityOfMouseOccupyingBin, 2)

        %%TO DO: from here 

        % conditional p of event given in bin is = (prob in bin + prob of event) / (prob in bin)

        % (conditional probability of observing a calcium events given mouse is in bin) *    
        % log ( (conditional probability of observing a calcium events given mouse is in bin) / (probability of observing k calcium events (0 or 1) for a cell
        jointProbability = neuronFiringProbability(1, i) .* probabilityOfMouseOccupyingBin(1, j);
        cProbEventGivenBin = jointProbability/probabilityOfMouseOccupyingBin(1, j); 
        
        SI_thisbin = cProbEventGivenBin*log(cProbEventGivenBin/neuronFiringProbability(1, i));
        
        MI_perbin(1, j) = SI_thisbin; 

    end

    MI_perCellperBin(i,:) = MI_perbin; 

MI_perCell = MI_perCellperBin * probabilityOfMouseOccupyingBin'; 
     
end
%% shuffle calcium traces and compute signal peaks 

%shuffle 
numShuffles = 100 ; 
cellTracesArray = table2array(cellTraces);
[numFrames, numNeurons] = size(cellTracesArray);
allPeaksShuffled = zeros(numFrames,numNeurons,numShuffles); 
for shuffle = 1:numShuffles
    disp(shuffle); 
    % Initialize an array to store the shuffled data
    shuffledCellTraces = cellTracesArray;
    % Loop through each neuron (column) and shuffle the frames
    for neuron = 1:numNeurons
        % Generate a random permutation of row indices
        shuffleIndices = randperm(numFrames);
        % Apply the shuffle to the current neuron (column)
        shuffledCellTraces(:, neuron) = cellTracesArray(shuffleIndices, neuron);
    end
    %
    [signalPeaksThisShuffle] = computeSignalPeaks(shuffledCellTraces, 'doMovAvg', 0, 'reportMidpoint', 1, 'numStdsForThresh', 2.5);
    allPeaksShuffled(:,:,shuffle)=signalPeaksThisShuffle; 
end 

%% calculate mutual information per cell

% overall firing probability for each neuron 
neuronFiringProbability = mean(signalPeaks, 1);

% probability of calcium event over occupany probability
expanded_probability_per_bin = repmat(probabilityOfMouseOccupyingBin, size(cellFiringProbabilityPerBin, 1), 1);
p_event_over_p_occupancy = cellFiringProbabilityPerBin ./ expanded_probability_per_bin;

% sum of event probabilities in each bin
MI_perCellperBin = zeros(size(cellFiringProbabilityPerBin,1), size(cellFiringProbabilityPerBin,2));
MI_perCell = zeros(1, size(cellFiringProbabilityPerBin,1));

% for each cell
size(cellFiringProbabilityPerBin,1)
for i = 1:size(cellFiringProbabilityPerBin,1)

    MI_perbin = zeros(1, size(probabilityOfMouseOccupyingBin, 2));
    
    %for each bin  
    for j = 1:size(probabilityOfMouseOccupyingBin, 2)

        %%TO DO: from here 

        % conditional p of event given in bin is = (prob in bin + prob of event) / (prob in bin)

        % (conditional probability of observing a calcium events given mouse is in bin) *    
        % log ( (conditional probability of observing a calcium events given mouse is in bin) / (probability of observing k calcium events (0 or 1) for a cell
        jointProbability = neuronFiringProbability(1, i) .* probabilityOfMouseOccupyingBin(1, j);
        cProbEventGivenBin = jointProbability/probabilityOfMouseOccupyingBin(1, j); 
        
        SI_thisbin = cProbEventGivenBin*log(cProbEventGivenBin/neuronFiringProbability(1, i));
        
        MI_perbin(1, j) = SI_thisbin; 

    end

    MI_perCellperBin(i,:) = MI_perbin; 

MI_perCell = MI_perCellperBin * probabilityOfMouseOccupyingBin'; 
     
end
























