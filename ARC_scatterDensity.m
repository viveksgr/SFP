function ARC_scatterDensity(x1, x2)
    % Check if the inputs are of equal length
    if length(x1) ~= length(x2)
        error('x1 and x2 must be vectors of the same length.');
    end

    % Calculate the point density
    [density, xi, yi] = ksdensity([x1, x2], [x1, x2]);
    
    % Create the scatter plot
    % figure;
    scatter(x1, x2, 15, density, 'filled'); % Adjust the size to 15 for visibility
    colormap(parula); % Use the jet colormap, or choose another if preferred
    colorbar; % Display a colorbar to indicate the density scale
    % title('Scatter Plot Colored by Density');
    % xlabel('x1');
    % ylabel('x2');
    axis tight; % Fit the axes tightly around the data
end
