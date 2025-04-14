function [p_accu, predictions] = ClassifyNestedPCA_SVM(nmatred, grp, nfolds)
% CLASSIFYNESTEDPCA_SVM Performs nested cross-validation with PCA feature selection
%                       and LIBSVM classification.
%
%   [p_accu, predictions] = CLASSIFYNESTEDPCA_SVM(nmatred, grp, nfolds)
%
%   INPUT:
%       nmatred : [N x P] data matrix, where N is number of trials,
%                 P is number of features.
%       grp     : [N x 1] integer or class labels for classification.
%       nfolds  : Number of folds in outer cross-validation (default = 10).
%
%   OUTPUT:
%       p_accu     : Overall performance metric (correlation or accuracy).
%       predictions : [N x 1] vector of predicted labels for each trial
%                     in a leave-one-out manner across outer folds.
%
%   LIBSVM is used for classification with RBF kernel by default. Adjust
%   parameters or kernel types as needed.
%
%   EXAMPLE:
%       % Suppose we have a data matrix "nmatred" of size [200 x 100]
%       % and a label vector "grp" of length 200. We run 10-fold outer CV:
%       [accu, preds] = ClassifyNestedPCA_SVM(nmatred, grp, 10);
%       disp(accu);

    if nargin < 3
        nfolds = 10;
    end

    % Outer cross-validation setup
    outerFolds = crossvalind('Kfold', length(grp), nfolds);

    predictions = zeros(size(grp));  % to store final test predictions

    % Range of possible PCs to try
    % Adjust as needed, e.g., 1:5:50 for big data, or a smaller range
    pcRange = 1:1:12;  % e.g. test up to 30 PCs in increments of 5

    for foldIdx = 1:nfolds
        % Outer train/test split
        testMask = (outerFolds == foldIdx);
        trainMask = ~testMask;

        trainX = nmatred(trainMask, :);
        trainY = grp(trainMask);
        testX  = nmatred(testMask, :);
        testY  = grp(testMask);

        % Inner cross-validation on train set for PCA dimension selection
        innerFolds = crossvalind('Kfold', length(trainY), nfolds);  % same folds or not necessarily

        bestPC = 1;
        bestPerf = 0;

        for npcs = pcRange
            perfVals = zeros(nfolds,1);

            for innerIdx = 1:nfolds
                innerTestMask = (innerFolds == innerIdx);
                innerTrainMask = ~innerTestMask;

                innerTrainX = trainX(innerTrainMask, :);
                innerTrainY = trainY(innerTrainMask);
                innerTestX  = trainX(innerTestMask, :);
                innerTestY  = trainY(innerTestMask);

                % PCA on innerTrainX
                [coeff, scoreTrain] = pca(innerTrainX);

                % Take first npcs components
                innerTrainPCs = scoreTrain(:, 1:npcs);

                % Transform innerTestX using same PCA components
                % 1) center
                mu = mean(innerTrainX,1);
                innerTestCentered = innerTestX - mu;

                % 2) project onto top PCs
                innerTestPCs = innerTestCentered * coeff(:, 1:npcs);

                % Train SVM on innerTrainPCs
                % Adjust kernel, cost as needed
                svmOpt = '-t 2 -c 1 -q'; % RBF kernel, cost=1
                model = svmtrain(innerTrainY, innerTrainPCs, svmOpt);

                % Predict
                preds = svmpredict(innerTestY, innerTestPCs, model, '-q');

                % Evaluate performance (accuracy here)
                acc = sum(preds == innerTestY) / length(innerTestY);
                perfVals(innerIdx) = acc;
            end

            meanPerf = mean(perfVals);
            if meanPerf > bestPerf
                bestPerf = meanPerf;
                bestPC = npcs;
            end
        end

        % We have bestPC for this fold. Retrain on entire fold's training set.
        [coeffFull, scoreFull] = pca(trainX);
        trainPCs = scoreFull(:, 1:bestPC);

        % Transform test set
        muTrain = mean(trainX,1);
        testCentered = testX - muTrain;
        testPCs = testCentered * coeffFull(:, 1:bestPC);

        % Train SVM on the entire training set
        svmOpt = '-t 2 -c 1 -q'; 
        finalModel = svmtrain(trainY, trainPCs, svmOpt);

        % Predict on the outer test set
        foldPreds = svmpredict(testY, testPCs, finalModel, '-q');
        predictions(testMask) = foldPreds;
    end

    % Overall performance metric
    % Could be accuracy or correlation
    % e.g. accuracy:
    p_accu = sum(predictions == grp) / length(grp);

    % or correlation, if numeric labels:
    % p_accu = fastcorr(predictions, grp);

end
