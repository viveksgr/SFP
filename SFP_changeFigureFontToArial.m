function SFP_changeFigureFontToArial(figHandle)
    % Ensure the figure handle is provided, if not use current figure
    if nargin < 1 || isempty(figHandle)
        figHandle = gcf; % Get current figure
    end

    % Find all text objects in the figure
    textObjects = findall(figHandle, 'Type', 'text');
    
    % Change font of all text objects to Arial
    set(textObjects, 'FontName', 'Arial');
    
    % Find all axes in the figure and change their font (this includes labels)
    axesObjects = findall(figHandle, 'Type', 'axes');
    set(axesObjects, 'FontName', 'Arial');
    
    % If the figure contains colorbars, change their font as well
    colorbarObjects = findall(figHandle, 'Type', 'colorbar');
    set(colorbarObjects, 'FontName', 'Arial');
    
    % Apply to title and labels explicitly in case they are missed
    for ax = axesObjects'
        title(ax, 'FontName', 'Arial');
        xlabel(ax, 'FontName', 'Arial');
        ylabel(ax, 'FontName', 'Arial');
        zlabel(ax, 'FontName', 'Arial'); % For 3D plots
    end
end
