function [wcss,idx,C] = SFP_kmeansElbowMethod(X, clusterRange,dister)
    if nargin<3
        dister = 'sqeuclidean';
    end
    
    if sum(isnan(X),'all')>0; warning('X has Nans. Replaced with zeros.'); end
    X(isnan(X))=0;
    % X is the Nxp data matrix with p features
    % clusterRange is a vector specifying the range of k values to test
    
    % Initialize the vector to hold the within-cluster sum of squares (WCSS) for each k
    wcss = zeros(length(clusterRange), 1);
    
    % Iterate over the range of cluster numbers
    for i = 1:length(clusterRange)
        k = clusterRange(i);
        
        % Perform k-means clustering with 'k' clusters
        [idx, C, sumd] = kmeans(X, k, 'MaxIter', 1000, 'Replicates', 1,'distance',dister);
        
        % Compute the total within-cluster sum of squares for this k
        wcss(i) = sum(sumd);
        
        % Display progress
        % fprintf('Completed k-means with k=%d clusters.\n', k);
    end
end