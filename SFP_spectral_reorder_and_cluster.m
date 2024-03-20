function [DO, i_sort, idx,sumd] = SFP_spectral_reorder_and_cluster(D, numClusters)
    if min(D(:)) < 0
        error('D must have positive elements only');
    end

    sz = size(D);
    Q = -D;
    Q(logical(eye(sz(1)))) = 0;
    Q(logical(eye(sz(1)))) = sum(-Q, 2);

    T = zeros(size(Q));
    T(logical(eye(sz(1)))) = 1./sqrt(sum(D,2));

    [V, ~] = eig(T*Q*T);
    [v, i_sort] = sort(V(:,2));

    DO = D(:,i_sort);
    DO = DO(i_sort,:);

    % Use k-means on the components of the second eigenvector for clustering
    % Applying k-means directly to the second eigenvector's components
    [idx,~,sumd] = kmeans(real(V(:,2)), numClusters, 'MaxIter', 1000, 'Replicates', 100);
end
