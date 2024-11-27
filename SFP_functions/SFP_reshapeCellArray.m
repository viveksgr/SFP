function result = SFP_reshapeCellArray(cellArray3D)

    % 
    % idx2 = [1:10 19 11:14 19 15 19 16:18]; % Use this is argsort in sub(2 and 3) to match the labels
    % idx1 = [1:18 19 19 19];

    % Extract the size of the input 3D cell array
    [n1, n2, n3] = size(cellArray3D);
    
    % Initialize the output cell array
    result = cell(n1, n2);
    
    % Loop over the first two dimensions
    for i = 1:n1
        for j = 1:n2
            % Initialize a list to hold the vectors for the current cell
            % vectorList = cell(1, n3-2);
             vectorList = cell(1, n3);
            
            % Extract vectors from the third dimension
            for k = 1:n3%-2
                % vectorList{k} = cellArray3D{i, j, k+2};
                vectorList{k} = cellArray3D{i, j};
            end
            
            % Concatenate the vectors horizontally to form a matrix
            result{i, j} = [vectorList{:}];
        end
    end
end
