function [alignedCells, indexCells] = SFP_GaleShapley(cellArray)

    % cellArray = permute(cellArray,[2 1 3]);

    n = length(cellArray);    % Number of matrices in the cell array
    firstMatrix = cellArray{1};
    
    p1 = size(firstMatrix, 1);
    p2 = size(firstMatrix, 2);
    
    % Initialize the aligned output and indices cells
    alignedCells = cell(n, 1);
    indexCells = cell(n-1, 1);
    
    % Keep the first matrix unchanged in the output
    alignedCells{1} = firstMatrix;
    
    % Loop through the remaining matrices in the cell array
    for i = 2:n
        currentMatrix = cellArray{i};
        % Initialize the preference matrices for Gale-Shapley
        men_pref = zeros(p2, p2);
        women_pref = zeros(p2, p2);
        
        % Calculate the correlation matrix and fill in preferences
        for j = 1:p2
            % Get the correlations of the j-th column of the first matrix with all columns of the current matrix
            correlations = max(corr(firstMatrix(:, j), currentMatrix),corr(-firstMatrix(:, j), currentMatrix));
            % Sort correlations to find the preference order (descending)
            [~, men_pref(j, :)] = sort(correlations, 'descend');
            % Similarly, each column of the current matrix prefers to align with the best match from the first matrix
            [~, sortedIdx] = sort(correlations, 'descend');
            women_pref(:, j) = sortedIdx;
        end
        
        % Use Gale-Shapley algorithm to find the stable matching
        stableMatch = galeshapley(p2, men_pref, women_pref);
        
        % Reorder columns of the current matrix according to the stable match and store them
        alignedCells{i} = currentMatrix(:, stableMatch);
        indexCells{i-1} = stableMatch;
    end

      % alignedCells = permute(cellArray,[2 1 3]);
       % indexCells = permute(cellArray,[2 1 3]);
end



function stablematch = galeshapley(N, men_pref, women_pref)

men_free = zeros(N,1);
women_suitor = zeros(N,N);
women_partner = zeros(N,1);
rank = zeros(N,N);


for i = 1:N
    for j = 1:N
        for k = 1:N
        if(women_pref(i,k) == j)
            rank(i,j) = k;
        end
        end
    end
end

while (min(women_partner) == 0)
    for i = 1:N
        if (men_free(i,1) == 0)
            next = find(men_pref(i,:) > 0, 1);
            women_suitor(men_pref(i,next),i) = i;
            men_pref(i,next) = 0;
        end
    end
    for i = 1:N
        for j = 1:N
            if(women_suitor(i,j) ~= 0)
                if(women_partner(i,1) == 0)
                    women_partner(i,1) = women_suitor(i,j);
                    men_free(j,1) = 1;
                end
                if(women_partner(i,1) ~= 0)
                if(rank(i,women_suitor(i,j)) < rank(i,women_partner(i,1)))
                    men_free(women_partner(i,1),1) = 0;
                    women_partner(i,1) = women_suitor(i,j);
                    men_free(j,1) = 1;
                    
                end
                end
            end
        end
    end
end

stablematch = women_partner;

end

          
        
        
        
        
