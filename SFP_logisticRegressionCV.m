function [avgAccuracy, pValue] = SFP_logisticRegressionCV(y, g)
    if ~iscategorical(y)
        y = categorical(y);
    end
    
    % Number of observations and categories
    N = numel(y);
    numCategories = numel(categories(y));
    
    % Initialize array to store test set predictions and actual labels
    allPredictions = zeros(N, 1);
    allTrueLabels = grp2idx(y); % Convert to numeric indices for comparison
    
    % Generate indices for 10-fold cross-validation
    cvIndices = crossvalind('Kfold', N, 10);
    
    % Perform 10-fold CV
    for k = 1:10
        testIdx = (cvIndices == k);
        trainIdx = ~testIdx;

        % Fit the multinomial logistic regression model on the training set
        B = mnrfit(g(trainIdx, :), y(trainIdx));
        
        % Predict on the test set
        yPredProb = mnrval(B, g(testIdx, :));
        [~, yPred] = max(yPredProb, [], 2);
        
        % Store predictions and true labels
        allPredictions(testIdx) = yPred;
    end
    
    % Compute overall accuracy
    avgAccuracy = mean(allPredictions == allTrueLabels);
    
    % Compute p-value using permutation test
    numPermutations = 1000;
    permutedAccuracies = zeros(numPermutations, 1);
    for i = 1:numPermutations
        permutedLabels = allTrueLabels(randperm(N));
        permutedAccuracies(i) = mean(allPredictions == permutedLabels);
    end
    
    % P-value: proportion of permuted accuracies >= observed accuracy
    pValue = mean(permutedAccuracies >= avgAccuracy);
    
    disp(['Average accuracy: ', num2str(avgAccuracy)]);
    disp(['P-value: ', num2str(pValue)]);
end
