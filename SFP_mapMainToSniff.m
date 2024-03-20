function indices = SFP_mapMainToSniff(onsets_main, onsets_sniff)
    % Initialize the output variable with zeros
    indices = zeros(size(onsets_main));
    
    % Iterate through each element in onsets_main
    for i = 1:length(onsets_main)
        % Check each cell in onsets_sniff
        for j = 1:length(onsets_sniff)
            % If the current element from onsets_main is found in the j-th cell of onsets_sniff
            if any(onsets_main(i) == onsets_sniff{j})
                % Assign the index of the cell to the corresponding position in indices
                indices(i) = j;
                % Break the loop once a match is found since onsets_main are distinct
                break;
            end
        end
    end
end
