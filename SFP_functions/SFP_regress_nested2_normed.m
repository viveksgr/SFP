function [p_accu, predictions, p_accut, bestC] = SFP_regress_nested2_normed(nmatred, grp, nfolds, nperm)
if nargin < 3
    nfolds = 10; % Default to 10-fold cross-validation if not specified
end
if nargin < 4
    nperm = 1; % Default to 1 permutation if not specified
end

nfeatures = size(nmatred,2);

accumat = zeros(nperm,1);
for pp = 1:nperm
    % Initialize variables
    predictions = zeros(size(grp));

    % Cross-validation setup
    outer_trainind = crossvalind('Kfold', length(grp), nfolds);

    for trl_idx = 1:nfolds
        % Define the test set for this fold
        testind_log = (outer_trainind == trl_idx);
        testX = nmatred(testind_log, :);
        testY = grp(testind_log);

        % Define the training set for this fold
        trainX = nmatred(~testind_log, :);
        trainY = grp(~testind_log);

        % Normalize training data
        mu = mean(trainX);
        sigma = std(trainX);
        trainX_normalized = (trainX - mu) ./ sigma;

        % Normalize test data using the same mean and std from training data
        testX_normalized = (testX - mu) ./ sigma;

        % Inner cross-validation for parameter tuning
        inner_trainind = crossvalind('Kfold', length(trainY), nfolds);
        bestC = 1; % Default C
        bestMSE = inf; % Track the best (lowest) MSE

        % Parameter grid for C (broader and more granular)
        C_vals = logspace(-2, 1, 7);  % Example: spans from 0.01 to 10 in logarithmic scale
        kernel_scales = 1/nfeatures*2 .^ (-6:2:6);
        % kernel_scales = 2 .^ (-6:2:6);

        for C = C_vals
            for kernel_scale = kernel_scales
                mse_inner = zeros(nfolds, 1);
                for inner_idx = 1:nfolds
                    inner_test_log = (inner_trainind == inner_idx);
                    inner_train_log = ~inner_test_log;

                    % Split inner training and validation sets
                    inner_trainX = trainX_normalized(inner_train_log, :);
                    inner_trainY = trainY(inner_train_log);
                    inner_testX = trainX_normalized(inner_test_log, :);
                    inner_testY = trainY(inner_test_log);

                    opt = sprintf('-t 2 -c %f -g %f -q', C,kernel_scale);
                    inner_mdl = svmtrain(inner_trainY, inner_trainX, opt);
                    inner_pred = svmpredict(inner_testY, inner_testX, inner_mdl, '-q');

                    % Calculate the cost function
                    % mse_inner(inner_idx) = 1 - fastcorr(inner_pred, inner_testY);
                    mse_inner(inner_idx) = 1 - sum(inner_pred==inner_testY)/length(inner_testY);
                end
                meanMSE = mean(mse_inner);
                if meanMSE < bestMSE
                    bestMSE = meanMSE;
                    bestC = C;
                    bestKernelScale = kernel_scale;
                end
            end
        end

        % Train model on the entire training set with best parameters found
        opt = sprintf('-t 2 -c %f -g %f -q', bestC,bestKernelScale);
        final_mdl = svmtrain(trainY, trainX_normalized, opt);
        final_predictions = svmpredict(testY, testX_normalized, final_mdl, '-q');

        predictions(testind_log) = final_predictions;
    end

    % Compute overall performance
    accumat(pp) = fastcorr(predictions, grp);
end

p_accu = mean(accumat);
[~, p_accut] = ARC_r2t(p_accu, length(grp));  % Assuming ARC_r2t is a function to compute p-value from correlation
end
