function [cosineSimilarity,mseLoss] = SFP_measureInformationLoss(originalData,reconstructedData)
    % Calculate variance across rows for each time point in the original data

     % Compute the Mean Squared Error (MSE) between original and reconstructed data
    mseLoss = mean((originalData - reconstructedData).^2, 1);
    % Display results
    % fprintf('Percentage of Variance Lost at Each Time Point:\n');
    % disp(infoLossPercent);

     % Normalize data to unit vectors
    originalNorm = sqrt(sum(originalData.^2, 1));
    reconstructedNorm = sqrt(sum(reconstructedData.^2, 1));
    % Compute cosine similarity
    cosineSimilarity = sum(originalData .* reconstructedData, 1) ./ (originalNorm .* reconstructedNorm);
end
