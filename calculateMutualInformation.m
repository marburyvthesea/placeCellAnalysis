function [MI_perCellThisShuffle, MI_perCellperBinThisShuffle] = calculateMutualInformation(cellFiringProbabilityPerBinThisShuffle, neuronFiringProbabilityThisShuffle, probabilityOfMouseOccupyingBin)
    % Calculate mutual information for each neuron and each spatial bin
    %
    % INPUTS:
    % cellFiringProbabilityPerBinThisShuffle - matrix of firing probabilities per bin (neurons x bins)
    % neuronFiringProbabilityThisShuffle - vector of overall firing probabilities per neuron (1 x neurons)
    % probabilityOfMouseOccupyingBin - vector of probabilities of mouse occupying each spatial bin (1 x bins)
    %
    % OUTPUTS:
    % MI_perCellThisShuffle - vector of mutual information per neuron (1 x neurons)
    % MI_perCellperBinThisShuffle - matrix of mutual information per neuron per bin (neurons x bins)

    % Initialize output matrices
    numNeurons = size(cellFiringProbabilityPerBinThisShuffle, 1);
    numBins = size(cellFiringProbabilityPerBinThisShuffle, 2);
    MI_perCellperBinThisShuffle = zeros(numNeurons, numBins);
    MI_perCellThisShuffle = zeros(1, numNeurons);

    % Loop over each neuron to calculate mutual information
    for i = 1:numNeurons
        MI_perbin = zeros(1, numBins);

        % For each bin
        for j = 1:numBins
            % Calculate joint probability (P(event and bin))
            jointProbability = neuronFiringProbabilityThisShuffle(i) * probabilityOfMouseOccupyingBin(j);
            
            % Conditional probability of event given bin (P(event | bin))
            if probabilityOfMouseOccupyingBin(j) > 0
                cProbEventGivenBin = jointProbability / probabilityOfMouseOccupyingBin(j);
            else
                cProbEventGivenBin = 0; % Avoid division by zero if bin occupancy is zero
            end

            % Calculate mutual information contribution for this bin
            if cProbEventGivenBin > 0 && neuronFiringProbabilityThisShuffle(i) > 0
                SI_thisbin = cProbEventGivenBin * log(cProbEventGivenBin / neuronFiringProbabilityThisShuffle(i));
            else
                SI_thisbin = 0; % Set to zero if any probability is zero to avoid log issues
            end

            MI_perbin(j) = SI_thisbin;
        end

        % Store the mutual information per bin for this neuron
        MI_perCellperBinThisShuffle(i, :) = MI_perbin;

        % Calculate and store the overall mutual information for this neuron
        MI_perCellThisShuffle(i) = sum(MI_perbin .* probabilityOfMouseOccupyingBin);
    end
end
