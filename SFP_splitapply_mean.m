function [S_omat_vals_r_reduced] = SFP_splitapply_mean(S_omat_vals_r,idx)
    % Number of clusters
    numClusters = length(unique(idx));

    % Step 2: Reduce S_omat_vals_r based on the clustering of Fless_mat_pruned
    [T,nvox] = size(S_omat_vals_r);
    S_omat_vals_r_reduced = zeros(numClusters,nvox);

    for i = 1:numClusters
        % Aggregate S_omat_vals_r based on cluster index
        % Note: Transposing idx to match the second dimension of S_omat_vals_r
        clusterMembers = S_omat_vals_r(idx == i,:);
        S_omat_vals_r_reduced(i,:) = mean(clusterMembers); % Mean across the trials in the cluster
    end
end
