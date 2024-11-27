function SFP_clearLargeVariables(threshold)
    % Set default threshold if not provided
    if nargin < 1
        threshold = 20; % Default threshold in MB
    end
    
    % Get a list of all variables in the caller's workspace
    vars = evalin('caller', 'who');
    
    % Iterate over each variable
    for i = 1:numel(vars)
        varName = vars{i};
        varSize = evalin('caller', sprintf('whos(''%s'').bytes', varName)) / (1024 * 1024); % Size in MB
        
        % Check if the variable size exceeds the threshold
        if varSize > threshold
            % Clear the variable from the caller's workspace
            evalin('caller', ['clear ', varName]);
            
            % Print a message indicating the removed variable and its size
            fprintf('Removed variable "%s" (%.2f MB)\n', varName, varSize);
        end
    end
end