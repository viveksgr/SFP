function VS_boxplot(cellArray)
    % Concatenate all cell data into a single vector
    dataVector = vertcat(cellArray{:});
    
    % Create a group vector indicating the group for each data point
    groupVector = arrayfun(@(i) repmat(i, size(cellArray{i}, 1), 1), ...
                           1:numel(cellArray), 'UniformOutput', false);
    groupVector = vertcat(groupVector{:});
    
    % Use boxplot to draw N boxes for the cell array
    boxplot(dataVector, groupVector, 'Notch', 'on', 'Labels', arrayfun(@num2str, 1:numel(cellArray), 'UniformOutput', false));
    
end
