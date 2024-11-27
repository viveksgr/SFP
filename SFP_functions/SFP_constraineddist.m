function [D2] = SFP_constraineddist(ZI, ZJ)
     % function [D2,sniffDists,perceptDists] = SFP_constraineddist(ZI, ZJ)
    % Partition the vectors into sniff and perceptual components
    alpha = 1;
    numSniffFeatures = 8; % Set this to the correct number of sniff features

    ZI_sniff = ZI(1:numSniffFeatures);
    ZI_percept = ZI(numSniffFeatures+1:end);

    ZJ_sniff = ZJ(:, 1:numSniffFeatures);
    ZJ_percept = ZJ(:, numSniffFeatures+1:end);

    % Calculate correlation-based "distances" for sniff data
    % % Assuming fastcorr outputs a correlation value between -1 and 1
    % sniffCorrs = arrayfun(@(row) fastcorr(ZI_sniff, ZJ_sniff(row, :)), 1:size(ZJ_sniff, 1));
    % 
    % % Convert correlations to distances
    % sniffDists = (1 - (sniffCorrs))/2+eps;

    % % Calculate Euclidean distances for sniff data
    sniffDists = sqrt(sum((ZJ_sniff - ZI_sniff).^2, 2))+eps;

    % Calculate correlation-based "distances" for perceptual data
    % perceptCorrs = arrayfun(@(row) fastcorr(ZI_percept, ZJ_percept(row, :)), 1:size(ZJ_percept, 1));
    % perceptCorrs = arrayfun(@(row) isequal(ZI_percept, ZJ_percept(row, :)), 1:size(ZJ_percept, 1));
    % 
    % % Convert correlations to distances
    % perceptDists = (1 - (perceptCorrs))/2;

    % Calculate Euclidean distances for perceptual data
    perceptDists = sqrt(sum((ZJ_percept - ZI_percept).^2, 2))+eps;

    perceptDists(perceptDists<0.001)=0.001;
    invperceptDists = 1 ./ perceptDists;
    perceptDistNonlinearAdjustment = alpha*tanh(invperceptDists*alpha);

    % Compute D2 as sniffDists + perceptDistNonlinearAdjustment
    D2 = sniffDists + perceptDistNonlinearAdjustment;


    % Combine distances: Here, simply averaging the distances. Adjust as necessary.
    % % D2 = regressmeout(sniffDists', perceptDists')';
    % D2 = ((1- alpha)*sniffDists -  alpha*perceptDists)+(alpha+eps);
end
