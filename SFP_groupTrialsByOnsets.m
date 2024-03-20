function [onsets,names,durations] = SFP_groupTrialsByOnsets(onsets_or, groupidx)
    % Find the unique group identifiers
    uniqueGroups = unique(groupidx);
    
    % Initialize the cell array to hold the onsets for each group
    onsets = cell(length(uniqueGroups), 1);
    names = cell(length(uniqueGroups), 1);
    durations = cell(length(uniqueGroups), 1);
    % Iterate through each unique group
    for i = 1:length(uniqueGroups)
        % Find indices of trials belonging to the current group
        currentGroupIdx = groupidx == uniqueGroups(i);
        
        % Extract the onsets for the current group and sort them
        onsets{i} = sort(onsets_or(currentGroupIdx));
        names{i} = num2str(i);
        durations{i} = 0;
    end
end
