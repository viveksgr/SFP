function [totalMI, pValue] = sfP_MI(MatrixA, MatrixB, numShuffles)
    [numClusters, numFeatures] = size(MatrixA);
    
    % Calculate original mutual information
    originalMI = zeros(1, numFeatures);
    for i = 1:numFeatures
        originalMI(i) = mutualInfoEstimate(MatrixA(:, i), MatrixB(:, i));
    end
    totalMI = sum(originalMI);  % Sum of MI across features
    
    % Perform shuffle test
    shuffleMIs = zeros(1, numShuffles);
    for n = 1:numShuffles
        shuffledB = MatrixB(randperm(numClusters), :);  % Shuffle rows of MatrixB
        shuffledMI = zeros(1, numFeatures);
        for i = 1:numFeatures
            shuffledMI(i) = mutualInfoEstimate(MatrixA(:, i), shuffledB(:, i));
        end
        shuffleMIs(n) = sum(shuffledMI);
    end
    
    % Calculate p-value
    pValue = mean(shuffleMIs >= totalMI);
    
    % Display results
    % fprintf('Original Total MI: %f\n', totalMI);
    % fprintf('P-value: %f\n', pValue);
end

function mi = mutualInfoEstimate(x, y)
    numBins = 20;
    jointHist = hist3([x(:), y(:)], {linspace(min(x), max(x), numBins), linspace(min(y), max(y), numBins)});
    jointProb = jointHist / sum(jointHist(:));
    marginalProbX = sum(jointProb, 2);
    marginalProbY = sum(jointProb, 1);
    
    hX = -sum(marginalProbX .* log(marginalProbX + eps));
    hY = -sum(marginalProbY .* log(marginalProbY + eps));
    hXY = -sum(jointProb(:) .* log(jointProb(:) + eps));
    
    mi = hX + hY - hXY;
end
