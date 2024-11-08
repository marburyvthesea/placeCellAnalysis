function cellFiringProbabilityPerBinThisShuffle = calculateFiringProbability(event_rateThisShuffle, x_position, bin_edges)
    % Calculate firing probability per bin for each neuron
    %
    % INPUTS:
    % event_rateThisShuffle - matrix of event rates (frames x neurons)
    % x_position - vector of x positions (frames x 1)
    % bin_edges - vector defining the edges of spatial bins
    %
    % OUTPUT:
    % cellFiringProbabilityPerBinThisShuffle - matrix of firing probabilities per bin (neurons x bins)

    % Transpose event_rateThisShuffle for ease of processing
    event_rate_t = event_rateThisShuffle';
    numNeurons = size(event_rate_t, 1);
    numBins = length(bin_edges) - 1; % Number of bins
    cellFiringProbabilityPerBinThisShuffle = zeros(numNeurons, numBins);

    % Loop over each neuron
    for neuron = 1:numNeurons
        % Find the indices where the event rate for this neuron is non-zero
        nonzero_indices = find(event_rate_t(neuron, :) ~= 0);
        
        % Get x positions at these indices
        x_position_when_active = x_position(nonzero_indices);
        
        % Calculate histogram counts for each bin (number of times the neuron is active in each bin)
        [event_counts_per_x_bin, ~] = histcounts(x_position_when_active, bin_edges);
        
        % Calculate the probability for each bin (fraction of active frames in each bin)
        cellFiringProbabilityPerBinThisShuffle(neuron, :) = event_counts_per_x_bin / length(event_rate_t(neuron, :));
    end
end