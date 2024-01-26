function processFiles(directory)
% Get all items in the directory
items = dir(directory);
names = {};
dateNums = [];
paths = {};
% Iterate through all items

for i = 1:length(items)
    item = items(i);

    % If item is a subdirectory and not '.' or '..'
    if item.isdir && ~strcmp(item.name, '.') && ~strcmp(item.name, '..') && ~strcmp(item.name, '.git')
        % Recursively process this subdirectory
        processFiles(fullfile(directory, item.name));
    else
        % If item is a MATLAB file
        if endsWith(item.name, '.m')
            % Extract file's date from its name
            parts = split(item.name, ' ');

            if length(parts) > 2
                dateStr = parts{2};
            else
                dateStr = '';
            end
            dateNum = item.datenum;
            % Convert the date string to a date number

            % Store the file's information
            names{end+1} = item.name;
            dateNums(end+1) = dateNum;
            paths{end+1} = fullfile(directory, item.name);    
        end
    end
end
filesInfo = struct('name', names, 'dateNum', num2cell(dateNums), 'path', paths);
    


% If there are files in this directory
if ~isempty(filesInfo)
    % Sort files based on their date numbers
    [~, indices] = sort([filesInfo.dateNum], 'descend');
    filesInfo = filesInfo(indices);

    % Rename files and keep only the most recent version
    previousFilenames = {};
    for i = 1:length(filesInfo)
        file = filesInfo(i);

        % Extract the base filename
        parts = split(file.name, ' ');
        baseFilename = parts{1};

        % If this file has a newer version, delete it
        if ismember(baseFilename, previousFilenames)
            delete(file.path);
        else
            % Rename the file to remove the date
            movefile(file.path, fullfile(directory, sprintf('%s.m',baseFilename)));
            previousFilenames{end+1} = baseFilename; %#ok<AGROW>
        end
    end
end
end

