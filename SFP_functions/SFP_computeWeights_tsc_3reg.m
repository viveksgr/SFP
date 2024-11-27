function [weights, t_scores] = SFP_computeWeights_tsc_3reg(x1, x2, x3, y,sw)

if nargin<4
    sw = false;
end

if sw
    % x1 = matrixToVecNoDiag(x1);
    % x2 = matrixToVecNoDiag(x2);
    % y = matrixToVecNoDiag(y);
    x1 = extractOffDiagonalBlock(x1);
    x2 = extractOffDiagonalBlock(x2);
    x3 = extractOffDiagonalBlock(x3);   
    y = extractOffDiagonalBlock(y);
end

[r,~] = find(isnan([x1,x2,x3,y]));
x1(r,:)=[];
x2(r,:)=[];
x3(r,:)=[];
y(r,:) = [];


% Combine the regressors into a design matrix, including a column for the intercept
X = [ones(length(x1), 1), x1(:), x2(:),x3(:)];

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
