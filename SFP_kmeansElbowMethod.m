function [wcss,idx,C] = SFP_kmeansElbowMethod(X, clusterRange)
    % X is the Nxp data matrix with p features
    % clusterRange is a vector specifying the range of k values to test
    
    % Initialize the vector to hold the within-cluster sum of squares (WCSS) for each k
    wcss = zeros(length(clusterRange), 1);
    
    % Iterate over the range of cluster numbers
    for i = 1:length(clusterRange)
        k = clusterRange(i);
        
        % Perform k-means clustering with 'k' clusters
        [idx, C, sumd] = kmeans(X, k, 'MaxIter', 1000, 'Replicates', 10);
        
        % Compute the total within-cluster sum of squares for this k
        wcss(i) = sum(sumd);
        
        % Display progress
        % fprintf('Completed k-means with k=%d clusters.\n', k);
    end
end