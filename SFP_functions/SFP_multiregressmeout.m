function [outputM] = SFP_multiregressmeout(M1, M2)
    % Ensure M1 and M2 have the same number of rows
    if size(M1, 1) ~= size(M2, 1)
        error('M1 and M2 must have the same number of rows.');
    end

    % Initialize the output matrix
    outputM = zeros(size(M1));

    % Iterate through each column of M1
    for colIdx = 1:size(M1, 2)
        % Select the column from M1
        vec1 = M1(:, colIdx);

        % Regress out each column of M2 from the selected column of M1
        % Create a design matrix for regression, including a constant term
        designMatrix = [ones(size(M2, 1), 1), M2];

        % Perform linear regression
        b = designMatrix \ vec1;

        % Calculate the residuals
        residuals = vec1 - designMatrix * b;

        % Store the residuals in the output matrix
        outputM(:, colIdx) = residuals;
    end
end
