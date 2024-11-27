function result = SFP_extractUpperTriangles(mainmat,utl_mask)
    % Get the size of the input matrix
    [N, P] = size(mainmat);
    
    % Calculate the number of elements in the upper triangle (excluding the diagonal)
    numUpperTriElements = N * (N - 1) / 2;
    
    % Initialize the result matrix with zeros
    result = zeros(numUpperTriElements, P);
    
    % Loop through each column of the mainmat matrix
    for col = 1:P
        % Compute the NxN matrix of absolute differences for the current column
        diffMatrix = -abs(mainmat(:, col) - mainmat(:, col)');
        
        % Convert the upper triangle to a column vector
        result(:, col) =  diffMatrix(utl_mask);
    end
end
