function [weights, t_scores] = ARC_multicomputeWeights_tsc(X, y, sw)

if nargin < 3
    sw = false;
end

if sw
    % Assuming extractOffDiagonalBlock function is designed to handle matrices where applicable
    X = extractOffDiagonalBlock(X);
    y = extractOffDiagonalBlock(y);
end

% Remove rows with NaN values in X or y
nanRows = any(isnan([X, y]), 2);
X(nanRows, :) = [];
y(nanRows) = [];

% % Standardize the input variables (X)
X = zscore(X);
X = [ones(size(X, 1), 1), X];

% % Standardize the output variable (y)
y = zscore(y);

% % Unnormalized:
% Solve for the weights
weights = X \ y;

% Calculate residuals
residuals = y - X * weights;

% Estimate the variance of the residuals
residual_variance = sum(residuals.^2) / (length(y) - size(X, 2));

% Calculate the standard error for each weight
standard_errors = sqrt(residual_variance * diag(inv(X' * X)));

% Calculate the t-scores
t_scores = weights ./ standard_errors;
end

function A2 = extractOffDiagonalBlock(A)
    % Ensure the input is a square matrix and its size is even
    [rows, cols] = size(A);
    if rows ~= cols || mod(rows, 2) ~= 0
        error('Input must be a square matrix with even size.');
    end

    % Calculate the size of the blocks
    blockSize = rows / 2;

    % Extract the off-diagonal block (e.g., top-right block)
    A2 = A(1:blockSize, (blockSize + 1):end);
    A2 = A2(:);

    % Alternatively, to extract the bottom-left block, use:
    % A2 = A((blockSize + 1):end, 1:blockSize);
end
