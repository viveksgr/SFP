function [idx_new, idx2_new] = sfp_decoupleClusteringIndices(idx, idx2, maxIterations, threshold)
    % Get the number of data points and unique cluster labels
    n = numel(idx);
    labels1 = unique(idx);
    labels2 = unique(idx2);
    
    % Initialize the new clustering indices
    idx_new = idx;
    idx2_new = idx2;
    
    % Initialize the iteration counter
    iter = 0;
    
    % Compute the initial adjusted Rand index
    ari = rand_index(idx_new, idx2_new, 'adjusted');
    
    % Iterate until the ARI is below the threshold or the maximum iterations are reached
    while ari > threshold && iter < maxIterations
        % Compute the co-occurrence matrix
        cooccurrence = sfp_cooccurrenceMatrix(idx_new, idx2_new);
        
        % Find the indices of the maximum value in the co-occurrence matrix
        [max_row, max_col] = find(cooccurrence == max(cooccurrence(:)), 1);
        
        % Find the trials corresponding to the maximum overlap
        overlap_mat = and((idx_new == labels1(max_row)),(idx2_new == labels2(max_col)));
        overlap_trials = find(overlap_mat);
        
        % Remove the overlapping trials from both clusterings
        idx_new(overlap_trials) = [];
        idx2_new(overlap_trials) = [];
        
        
        % Update the cluster labels
        labels1 = unique(idx_new);
        labels2 = unique(idx2_new);
        
        % Update the adjusted Rand index
        ari = rand_index(idx_new, idx2_new, 'adjusted');
        
        % Increment the iteration counter
        iter = iter + 1;
        if iter == maxIterations
            fprintf('Max limit reached')
            break
        end
    end
end

