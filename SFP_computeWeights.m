function weights = SFP_computeWeights(x1, x2, y)
    % Combine the regressors into a design matrix, including a column for the intercept
    X = [ones(length(x1), 1), x1(:), x2(:)];
    
    % Solve for the weights
    weights = X \ y(:);
end
