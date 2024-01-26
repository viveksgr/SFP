function [weights, t_scores] = SFP_computeWeights_tsc_4reg(x1, x2, x3, x4, y)

    x1 = matrixToVecNoDiag(x1);
    x2 = matrixToVecNoDiag(x2);
    x3 = matrixToVecNoDiag(x3);
    x4 = matrixToVecNoDiag(x4);
    y = matrixToVecNoDiag(y);    


    % Remove rows with NaN values in any of the variables
    [r,~] = find(isnan([x1, x2, x3, x4, y]));
    x1(r,:) = [];
    x2(r,:) = [];
    x3(r,:) = [];
    x4(r,:) = [];
    y(r,:)  = [];

    % Combine the regressors into a design matrix, including a column for the intercept
    X = [ones(length(x1), 1), x1(:), x2(:), x3(:), x4(:)];
    
    % Solve for the weights
    weights = X \ y(:);

    % Calculate residuals
    residuals = y(:) - X * weights;

    % Estimate the variance of the residuals
    residual_variance = sum(residuals.^2) / (length(y) - length(weights));

    % Calculate the standard error for each weight
    standard_errors = sqrt(residual_variance * diag(inv(X'*X)));

    % Calculate the t-scores
    t_scores = weights ./ standard_errors;
end


function vec = matrixToVecNoDiag(A)
    % Ensure the input is a square matrix
    [rows, cols] = size(A);
    if rows ~= cols
        error('Input must be a square matrix.');
    end

    % Convert the matrix to a vector
    vec = A(:);

    % Create an index for the diagonal elements
    diagIndex = 1:(rows+1):length(vec);

    % Remove the diagonal elements
    vec(diagIndex) = [];
end