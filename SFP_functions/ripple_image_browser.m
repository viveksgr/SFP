function ripple_image_browser(filepath, fileextension)
    % filepath: string of filepath with images: e.g. 'C:\Work\...'
    % fileextension: extension of images, e.g. jpg, or png

    % Get list of all images in the specified folder with the given extension
    imageFiles = dir(fullfile(filepath, strcat('*.', fileextension)));
    numImages = length(imageFiles);
    if numImages == 0
        error('No images found with the specified extension in the specified directory.');
    end

    % Initialize the figure
    hFig = figure('Name', 'Image Browser', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
    hAx = axes('Parent', hFig, 'Position', [0.1, 0.2, 0.8, 0.7]);
    hImage = imshow(fullfile(imageFiles(1).folder, imageFiles(1).name), 'Parent', hAx, 'InitialMagnification', 'fit');

    % Initialize index and response structure
    idx = 1; % Index of the current image being displayed
    responses = struct('FileName', {imageFiles.name}, 'ButtonPress', []); % Structure to store the file names and button responses

     % Create buttons
    hPrevButton = uicontrol('Style', 'pushbutton', 'String', 'Previous', 'Position', [150, 20, 100, 40], 'Callback', @(~,~) showImage(-1));
    hNextButton = uicontrol('Style', 'pushbutton', 'String', 'Next', 'Position', [550, 20, 100, 40], 'Callback', @(~,~) showImage(1));
    uicontrol('Style', 'pushbutton', 'String', 'Ripple', 'Position', [270, 20, 80, 40], 'Callback', @(~,~) recordAndMove('Ripple'));
    uicontrol('Style', 'pushbutton', 'String', 'HFO', 'Position', [360, 20, 80, 40], 'Callback', @(~,~) recordAndMove('HFO'));
    uicontrol('Style', 'pushbutton', 'String', 'Unsure', 'Position', [450, 20, 80, 40], 'Callback', @(~,~) recordAndMove('Unsure'));
    uicontrol('Style', 'pushbutton', 'String', 'Save', 'Position', [690, 20, 80, 40], 'Callback', @saveResponses);

    % Initialize button states
    updateButtonStates();

    function showImage(step)
        idx = idx + step;
        if idx < 1
            idx = 1;
        elseif idx > numImages
            idx = numImages;
        end
        set(hImage, 'CData', imread(fullfile(imageFiles(idx).folder, imageFiles(idx).name)));
        updateButtonStates();
    end

    function updateButtonStates()
        if idx <= 1
            set(hPrevButton, 'Enable', 'off');
        else
            set(hPrevButton, 'Enable', 'on');
        end
        
        if idx >= numImages
            set(hNextButton, 'Enable', 'off');
        else
            set(hNextButton, 'Enable', 'on');
        end
    end

    function recordAndMove(response)
        responses(idx).ButtonPress = response;
        if idx < numImages
            showImage(1);
        end
    end

    function saveResponses(~,~)
        save('responses.mat', 'responses');
        msgbox('Responses saved successfully!');
    end
end
