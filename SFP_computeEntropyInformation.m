function information = SFP_computeEntropyInformation(data)
 % Normalize the data to form a probability distribution for each column
    probMatrix = max(data, 0); % Ensure no negative values
    sumProb = sum(probMatrix, 1);
    probMatrix = bsxfun(@rdivide, probMatrix, sumProb); % Normalize columns to sum to 1
    
    % Define entropy function considering zero probabilities
    entropyFunc = @(p) -nansum(p(p > 0) .* log2(p(p > 0) + eps));
    
    % Compute entropy for each time point
    information = arrayfun(@(t) entropyFunc(probMatrix(:, t)), 1:size(data, 2));
end
