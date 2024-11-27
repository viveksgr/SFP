function SFP_colorfulLinePlot(dataVector)
    % Create a figure
    figure;

    % Get the number of data points
    numPoints = length(dataVector);

    % Generate x values based on the vector length
    xValues = 1:numPoints;

    % Get the 'parula' colormap with as many colors as there are points
    colors = parula(numPoints);

    % Hold on to allow multiple plots in the same figure
    hold on;

    % Plot each point individually
    for i = 1:numPoints
        plot(xValues(i:i+1), dataVector(i:i+1), 'Color', colors(i,:), 'LineWidth', 2);
    end

    % Hold off after plotting
    hold off;

    % Enhance the plot
    title('Colorful Line Plot');
    xlabel('Index');
    ylabel('Value');
    grid on;

    % Adjust axis for better visibility
    axis tight;
end
