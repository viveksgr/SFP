function ARC_scatterDensity(x1, x2,sw)


    if nargin<3
        sw = true; % Plot linear regression
    end
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

    
    hold on
    wt = [ones(size(x1)) x1]\x2;
    ypred = [ones(size(x1)) x1]*wt;
    plot(x1,ypred,'k')
    [~,tsc] = ARC_multicomputeWeights_tsc(x1, x2); 
    title(sprintf('m: %.3f; r: %.3f, p: %.3f',wt(2),fastcorr(x1,x2),1-tcdf(abs(tsc(2)),length(x1))))
end
