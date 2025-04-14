function [p_accu,predictions] = Classify_Permute_VS2( nmatred, grp, nfolds)


% SVM decoding 
% 
% Input
%   nmatred, trial x feature
%   grp, trial x1 
%   nfolds = folds of crossvalidation (default 10)
% 
% Output
%   p_accu, permutated accuracy

if nargin<3
    nfolds = 10;
end

% assert(length(grp)==size(nmatred,1),'Size of labels does not match training data')
trainind = crossvalind('Kfold',length(grp),nfolds);

pred_val = zeros(nfolds,1);
predictions = zeros(size(grp));
for trl_idx = 1 : nfolds
    testind_log = trainind==trl_idx;

    trainX = nmatred(~testind_log,:);
    testX = nmatred(testind_log,:);
    trainY = grp(~testind_log);
    testY = grp(testind_log);

    % Classification
    mdl = svmtrain( trainY,  trainX, ' -t 2 -c 1 -q');

    % Regression:
    % mdl = svmtrain( trainY,  trainX, '-s 3 -t 2 -c 1 -q');
    % eval( ['mdl = svmtrain( l1( :, 1), mat( l1( :, 2), :), '' ', opt.svm, ''');']);
    tmp = svmpredict(  testY , testX ,  mdl, ' -q ');
    predictions(testind_log) = tmp;
    pred_val(trl_idx) = sum(tmp==testY)/length(testY);
end % leave-one-out
% p_accu=mean(pred_val);
p_accu = fastcorr(predictions,grp);
% p_accu = sum(predictions == grp) / length(grp);

