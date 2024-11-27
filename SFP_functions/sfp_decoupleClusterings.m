function [idx1_new, idx2_new] = sfp_decoupleClusterings(idx1, idx2, tolerance)
    % Get the number of data points and unique cluster labels
    n = length(idx1);
    labels1 = unique(idx1);
    labels2 = unique(idx2);
    
    % Initialize the new clusterings
    idx1_new = idx1;
    idx2_new = idx2;
    
    % Iterate until the ARI between idx1_new and idx2_new is less than the
    % tolerance
    max_counter = 10000; 
    mp = 0;
    while (rand_index(idx1_new, idx2_new, 'adjusted') >= tolerance)
        mp = mp+1;
        % Randomly select a data point
        i = randi(n);
        
        % Find the current cluster labels of the selected data point
        label1 = idx1_new(i);
        label2 = idx2_new(i);
        
        % Find the data points in the same cluster as the selected data point
        mask1 = (idx1_new == label1);
        mask2 = (idx2_new == label2);
        
        % Randomly select a new cluster label for idx1_new
        new_label1 = labels1(randi(length(labels1)));
        while new_label1 == label1
            new_label1 = labels1(randi(length(labels1)));
        end
        
        % Assign the new cluster label to the data points in the same cluster as the selected data point
        idx1_new(mask1) = new_label1;
        
        % Randomly select a new cluster label for idx2_new
        new_label2 = labels2(randi(length(labels2)));
        while new_label2 == label2
            new_label2 = labels2(randi(length(labels2)));
        end
        
        % Assign the new cluster label to the data points in the same cluster as the selected data point
        idx2_new(mask2) = new_label2;

        if mp==max_counter-1
            fprintf('Max Iterations Reached...')
            break
        end
    end
end