function fileMatrix = SFP_organizeFiles(folderPath)
    % Get all .nii files in the folder
    files = dir(fullfile(folderPath, 'wSFP*.nii'));
    filenames = {files.name};

    % Initialize containers for numbers and strings
    numbers = zeros(length(filenames), 1);
    strings = cell(length(filenames), 1);

    % Parse filenames to extract numbers and strings
    for i = 1:length(filenames)
        tokens = regexp(filenames{i}, 'wSFP(\d+)_(.*).nii', 'tokens');
        numbers(i) = str2double(tokens{1}{1});  % Convert number part to double
        strings{i} = tokens{1}{2};  % Extract string part
    end

    % Find unique strings and sort them
    uniqueStrings = unique(strings);
    numStrings = length(uniqueStrings);

    % Create a cell array to hold sorted filenames
    fileMatrix = cell(3, numStrings);

    % Fill the cell array with filenames sorted by number and string
    for i = 1:length(filenames)
        % Find the index for each unique string
        stringIndex = find(strcmp(strings{i}, uniqueStrings));
        
        % Place filename in the correct spot in the cell array
        fileMatrix{numbers(i), stringIndex} = filenames{i};
    end
end
