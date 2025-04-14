function condTable = sfp_buildConditionTable(list_id, conditionIDs, condMap)
% BUILDCONDITIONTABLE Creates a table with group, condition-ID, and condition-name.
%
%   condTable = buildConditionTable(list_id, conditionIDs, condMap)
%
%   INPUTS:
%       list_id      : A 10x20 matrix of serial indices (1..160). 
%                     Each of the 20 columns is a group of 10 conditions.
%       conditionIDs : A vector of length 160 mapping serial indices 
%                     (1..160) -> numeric keys (e.g. 101, 105, etc).
%       condMap      : A containers.Map from numeric keys -> condition names.
%
%   OUTPUT:
%       condTable    : A table with columns (Group, ConditionID, ConditionName),
%                     stacked for all (10 x 20) entries.

    % Dimensions based on list_id
    [numConds, numGroups] = size(list_id);  % e.g., 10 rows, 20 columns
    totalRows = numConds * numGroups;       % e.g., 10 * 20 = 200

    % Preallocate arrays for storing table data
    groupArray      = zeros(totalRows, 1);
    conditionIDArray = zeros(totalRows, 1);
    conditionNameArray = strings(totalRows, 1);  % Use a string array

    rowCounter = 1;

    % Loop over each group (column) and each row (condition) in list_id
    for g = 1:numGroups
        for r = 1:numConds
            % 1. Get the serial index from list_id
            serialIndex = list_id(r, g);  % e.g. a value in [1..160]

            % 2. Map that serial index to the numeric key
            %    (assuming conditionIDs is a 1D array of length 160)
            numericKey = conditionIDs(serialIndex); 

            % 3. Retrieve the name from the map container
            condName = condMap(numericKey);

            % 4. Store data in the arrays
            groupArray(rowCounter) = g;
            conditionIDArray(rowCounter) = numericKey;
            conditionNameArray(rowCounter) = condName;

            rowCounter = rowCounter + 1;
        end
    end

    % Build a table from the arrays
    condTable = table(groupArray, conditionIDArray, conditionNameArray, ...
        'VariableNames', {'Group', 'ConditionID', 'ConditionName'});
end
