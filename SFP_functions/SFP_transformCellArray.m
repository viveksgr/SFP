function [diffMat1, diffMat2, diffMat3] = SFP_transformCellArray(dataCell)
    % Initialize cell arrays for the differences
    diffMat1 = cell(size(dataCell, 1), size(dataCell, 2));  % Q1 - Q3
    diffMat2 = cell(size(dataCell, 1), size(dataCell, 2));  % Q2 - Q3
    diffMat3 = cell(size(dataCell, 1), size(dataCell, 2));  % Q1 - Q2

    % Iterate over each cell in the dataCell array
    for i = 1:size(dataCell, 1)  % Loop over the first dimension (subjects)
        for j = 1:size(dataCell, 2)  % Loop over the second dimension (anatomical areas)
            if ~isempty(dataCell{i, j, 1}) && ~isempty(dataCell{i, j, 3})
                diffMat1{i, j} = dataCell{i, j, 1} - dataCell{i, j, 3};  % Q1 - Q3
            end
            if ~isempty(dataCell{i, j, 2}) && ~isempty(dataCell{i, j, 3})
                diffMat2{i, j} = dataCell{i, j, 2} - dataCell{i, j, 3};  % Q2 - Q3
            end
            if ~isempty(dataCell{i, j, 1}) && ~isempty(dataCell{i, j, 2})
                diffMat3{i, j} = dataCell{i, j, 1} - dataCell{i, j, 2};  % Q1 - Q2
            end
        end
    end
end
