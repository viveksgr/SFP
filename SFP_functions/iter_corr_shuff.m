function [corrmat] = iter_corr_shuff(pred_vox,test_vox)
nperm = 1000;
corrmat = zeros(size(pred_vox,1),nperm);
for tt = 1:1000
    shuff_vec = randperm(size(pred_vox,1));
    pred_vox = pred_vox(shuff_vec,:);

    % Make each column zero-mean
    pred_vox = bsxfun(@minus,pred_vox,mean(pred_vox,2));
    test_vox = bsxfun(@minus,test_vox,mean(test_vox,2));

    % L2 normalize each column
    pred_vox = bsxfun(@times,pred_vox,1./sqrt(sum(pred_vox.^2,2)));
    test_vox = bsxfun(@times,test_vox,1./sqrt(sum(test_vox.^2,2)));

    % Take the dot product of the columns and then sum
    corrmat(:,tt)=sum(pred_vox.*test_vox,2);
end

end

