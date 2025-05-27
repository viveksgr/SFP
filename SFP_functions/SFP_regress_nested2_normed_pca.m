function [p_accu, predictions, p_accut, bestC] = ...
         SFP_regress_nested2_normed_pca(nmatred, grp, nfolds, nperm)

if nargin < 3,  nfolds = 10;  end           % outer folds
if nargin < 4,  nperm  = 1;   end           % permutations

accumat = zeros(nperm,1);

for pp = 1:nperm
    
    predictions = zeros(size(grp));               % final test-set predictions
    outerFoldId = crossvalind('Kfold', length(grp), nfolds);

    for fold = 1:nfolds
        % ---------- split outer train / test ----------
        testMask  =  outerFoldId == fold;
        trainMask = ~testMask;

        trainX = nmatred(trainMask,:);   trainY = grp(trainMask);
        testX  = nmatred(testMask ,:);   testY  = grp(testMask);

        % ---------- z-score (row-wise normalisation) ----------
        mu    = mean(trainX,1);
        sigma = std(trainX,[],1);
        trainXn = (trainX - mu) ./ sigma;
        testXn  = (testX  - mu) ./ sigma;

        % ▶ ---------- PCA ON TRAINING DATA  (70 % VAR) ----------
        [coeff,score,latent,~,explained] = pca(trainXn,'Centered',false);
        cumVar = cumsum(explained);
        k      = find(cumVar >= 60, 1);               % #PCs to keep (≥70 %)
        trainPC = score(:,1:k);                       % training in PCA space
        testPC  = testXn * coeff(:,1:k);              % project test with same PCs
        nFeatPCA = k;                                 % dimensionality after PCA
        % --------------------------------------------------------

        % ---------- inner CV for SVM hyper-parameters ------------
        innerFoldId = crossvalind('Kfold', length(trainY), nfolds);
        C_vals       = logspace(-2,1,7);
        kernel_scales_base = 1/nFeatPCA;              % make γ inversely ∝ dim
        gamma_vals   = kernel_scales_base * 2.^(-6:2:6);

        bestC = C_vals(1);  bestG = gamma_vals(1);  bestErr = inf;

        for C = C_vals
            for G = gamma_vals
                errs = zeros(nfolds,1);
                for iF = 1:nfolds
                    valMask  = innerFoldId == iF;
                    trMask   = ~valMask;

                    trX = trainPC(trMask ,:); trY = trainY(trMask);
                    vlX = trainPC(valMask,:); vlY = trainY(valMask);

                    model = svmtrain(trY, trX, ...
                            sprintf('-t 2 -c %g -g %g -q', C, G));
                    pred  = svmpredict(vlY, vlX, model, '-q');
                    errs(iF) = 1 - sum(pred==vlY)/numel(vlY);  % classification error
                end
                mErr = mean(errs);
                if mErr < bestErr
                    bestErr = mErr;  bestC = C;  bestG = G;
                end
            end
        end

        % ---------- retrain on FULL outer-train & test -----------
        finalModel = svmtrain(trainY, trainPC, ...
                     sprintf('-t 2 -c %g -g %g -q', bestC, bestG));
        preds = svmpredict(testY, testPC, finalModel, '-q');
        predictions(testMask) = preds;
    end

    % one permutation’s outer-fold performance (accuracy here)
    accumat(pp) = sum(predictions==grp)/length(grp);
end

p_accu = mean(accumat);                        % grand accuracy across permutations
[~,p_accut] = ARC_r2t(p_accu,length(grp));    % user-supplied t-test helper
end
