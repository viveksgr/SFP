function [avgAccuracy, pValue,featureRMS] = SFP_logisticRegressionCV(y, g)
    if ~iscategorical(y)
        y = categorical(y);
    end
    
    % Number of observations and categories
    N = numel(y);
    numCategories = numel(categories(y));
    numFeatures = size(g,2);
    
    % Initialize array to store test set predictions and actual labels
    allPredictions = zeros(N, 1);
    allTrueLabels = grp2idx(y); % Convert to numeric indices for comparison
    
    % Generate indices for 10-fold cross-validation
    cvIndices = crossvalind('Kfold', N, 10);
    
    % Perform 10-fold CV
    nfolds = 10;
    betaw = zeros((numCategories-1)*(size(g,2)+1),nfolds);
    for k = 1:nfolds
        testIdx = (cvIndices == k);
        trainIdx = ~testIdx;
        
        X = g(trainIdx, :);
        y_vec = y(trainIdx);
        % Fit the multinomial logistic regression model on the training set
        B = fitmnr(X, y_vec);
        
        betaw(:,k) = table2array(B.Coefficients(:,"Value"));
        % Predict on the test set

        yPredProb = predict(B,g(testIdx, :));
        
        % Store predictions and true labels
        allPredictions(testIdx) = yPredProb;
    end
        
    % Initialize vector to store RMS values
    featureRMS = zeros(numFeatures, 1);

    % Compute RMS for each feature across all classes and folds
    for i = 1:numFeatures
        % Extract the weights for the i-th feature across all classes and folds
        % Skipping the first row if it includes intercepts, adjust indexing if different
        featureWeights = betaw(i+1:numFeatures+1:end, :); % Adjust if intercepts are at a different position
        featureRMS(i) = sqrt(mean(featureWeights(:).^2));
    end

    % Compute overall accuracy
    avgAccuracy = mean(allPredictions == allTrueLabels);
    
    % Compute p-value using permutation test
    % numPermutations = 1000;
    % permutedAccuracies = zeros(numPermutations, 1);
    % for i = 1:numPermutations
    %     permutedLabels = allTrueLabels(randperm(N));
    %     permutedAccuracies(i) = mean(allPredictions == permutedLabels);
    % end
    % 
    % % P-value: proportion of permuted accuracies >= observed accuracy
    % pValue = mean(permutedAccuracies >= avgAccuracy);

    N = numel(y);
    numClasses = numel(categories(y));
    chanceAccuracy = 1 / numClasses;

    % Calculate the number of correct predictions
    numCorrect = sum(allPredictions == allTrueLabels);
    % Calculate the p-value using the binomial distribution
    pValue = 1 - binocdf(numCorrect - 1, N, chanceAccuracy);
    
    disp(['Average accuracy: ', num2str(avgAccuracy)]);
    disp(['P-value: ', num2str(pValue)]);
end
