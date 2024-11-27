function [idx_new, idx2_new,ari] = sfp_decoupleClusteringIndices_prob(idx, idx2, A1_corr, A2_corr, maxIterations, threshold,kernel_k)

if nargin<7
    kernel_k = floor(length(unique(idx))/6);
end

    n = numel(idx);
    labels1 = unique(idx);
    labels2 = unique(idx2);
    
    idx_new = idx;
    idx2_new = idx2;
    
    iter = 0;
    ari = rand_index(idx_new, idx2_new, 'adjusted');
    
    while ari > threshold && iter < maxIterations
        cooccurrence = sfp_cooccurrenceMatrix(idx_new, idx2_new);
        
        max_cooc = max(cooccurrence(:));
        if max_cooc == 1
            fprintf('Maximum co-occurrence is 1, breaking loop.\n');
            break;  % Exit the loop if the maximum co-occurrence is 1
        end
        
        mean_cooc = mean(cooccurrence(:));
        
        [max_row, max_col] = find(cooccurrence == max_cooc, 1);
        overlap_mat = and((idx_new == labels1(max_row)), (idx2_new == labels2(max_col)));
        overlap_trials = find(overlap_mat);
        
        % Use a proportion of the max cooccurrence, aiming to reduce it but leave some overlap
        num_trials_to_reassign = floor(max_cooc - mean_cooc);
        num_trials_to_reassign = min(num_trials_to_reassign, length(overlap_trials));  % Ensure not to exceed available trials
        
        trials_to_reassign = randsample(overlap_trials, num_trials_to_reassign);

        % Reassign these trials based on similarity matrices
        current_cluster_idx1 = find(labels1 == idx_new(trials_to_reassign(1)));
        similarities1 = A1_corr(current_cluster_idx1, :);
        similarities1 = (similarities1+1)/2;
        similarities1(current_cluster_idx1) = 0;  % Remove the current cluster from possibilities
        probabilities1 = similarities1 / sum(similarities1);  % Normalize to create a probability distribution
        % Select top-k clusters
        
        [top_probs1, top_indices1] = maxk(probabilities1, kernel_k);
        top_probs1 = top_probs1 / sum(top_probs1);  % Normalize top-k probabilities


        current_cluster_idx2 = find(labels2 == idx2_new(trials_to_reassign(1)));
        similarities2 = A2_corr(current_cluster_idx2, :);
        similarities2 = (similarities2+1)/2;
        similarities2(current_cluster_idx2) = 0;  % Remove the current cluster from possibilities
        probabilities2 = similarities2 / sum(similarities2);  % Normalize to create a probability distribution
        [top_probs2, top_indices2] = maxk(probabilities2, kernel_k);
        top_probs2 = top_probs2 / sum(top_probs2);  % Normalize top-k probabilities

     
        for i = trials_to_reassign'
            new_cluster_idx1 = randsample(labels1(top_indices1), 1, true, top_probs1);
            idx_new(i) = new_cluster_idx1;

            new_cluster_idx2 = randsample(labels2(top_indices2), 1, true, top_probs2);
            idx2_new(i) = new_cluster_idx2;
        end
        
        ari = rand_index(idx_new, idx2_new, 'adjusted');
        iter = iter + 1;
        
        if iter == maxIterations
            fprintf('Max limit reached\n');
            break;
        end
    end
    ari(2)=rand_index(idx_new, idx, 'adjusted');
    ari(3)=rand_index(idx2_new, idx2, 'adjusted');
end
