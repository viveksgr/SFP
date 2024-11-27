function SFP_plotSortedLines(data)
    nr = size(data,1);
    nc = size(data,2);

    % Find the indices of the maxima in each row
    [~, maxIndices] = max(data, [], 2);
    
    % Sort the rows based on the location of the maxima
    [~, sortOrder] = sort(maxIndices);
    sortedData = data(sortOrder, :);
    
    % Create the figure
    figure;
    hold on; % Hold on to plot all lines on the same graph
    
    % Get a color map with 10 colors
    colors = parula( nr); % 'jet' is a nice gradient; other options could be 'hot', 'cool', etc.
    
    % Plot each row with a color corresponding to its sorted order
    for i = 1: nr
        plot(linspace(0,4,nc),sortedData(i, :), 'Color', colors(i, :), 'LineWidth', 2);
    end
    
    % Enhancements for better visualization
    title('Sniff traces: sniff clustering');
    xlabel('Time(s)');
    ylabel('Sniff flow rate');
    axis tight; % Fit the axes tightly to the data
    grid on; % Turn on the grid for better readability
    hold off; % Release the hold on current figure
end
