function [reducedData, V, explainedVariance] = SFP_temporalPCA(dataMatrix, varianceThreshold)
    if nargin<2
        varianceThreshold = 0.95;
    end


    % Input dataMatrix is N x T
    % Center the data by subtracting the mean of each time point
    centeredData = dataMatrix - mean(dataMatrix, 1);

    % Calculate the covariance matrix of the transposed data
    covMatrix = cov(centeredData);

    % Perform eigenvalue decomposition
    [V, D] = eig(covMatrix);
    
    % Sort eigenvalues and eigenvectors
    [d, ind] = sort(diag(D), 'descend');
    V = V(:, ind);
    D = diag(d);

    % Calculate explained variance
    totalVariance = sum(diag(D));
    explainedVariance = cumsum(diag(D)) / totalVariance;
    
    % Determine the number of components needed to explain the desired threshold of variance
    numComponents = find(explainedVariance >= varianceThreshold, 1, 'first');

    % Project the data onto the principal components
    reducedData = centeredData * V(:, 1:numComponents);
    % V = V(:, 1:numComponents);
    % D = D(1:numComponents, 1:numComponents);

    % Optionally return explained variance per component
    % explainedVariance = explainedVariance(1:numComponents);
end
