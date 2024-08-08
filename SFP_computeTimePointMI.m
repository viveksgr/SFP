function miValues = SFP_computeTimePointMI(originalData, reconstructedData, numBins)
    % Number of time points
    T = size(originalData, 2);
    miValues = zeros(1, T);
    
    % Process each time point separately
    for t = 1:T
        % Get the data for this time point
        orig = originalData(:, t);
        recon = reconstructedData(:, t);
        
        % Bin the data
        edges = linspace(min([orig; recon]), max([orig; recon]), numBins + 1);
        orig_binned = discretize(orig, edges);
        recon_binned = discretize(recon, edges);
        
        % Compute joint histogram
        jointHist = accumarray([orig_binned recon_binned], 1, [numBins numBins]);
        jointProb = jointHist / sum(jointHist(:));
        
        % Marginal probabilities
        probOrig = sum(jointProb, 2);
        probRecon = sum(jointProb, 1);
        
        % Calculate mutual information
        MI = 0;
        for i = 1:numBins
            for j = 1:numBins
                if jointProb(i,j) > 0  % Avoid log(0)
                    MI = MI + jointProb(i,j) * log2(jointProb(i,j) / (probOrig(i) * probRecon(j)));
                end
            end
        end
        
        % Store the MI for this time point
        miValues(t) = MI;
    end
end
