% function maxCorr = correlateDistanceMatrix(N, P)

    p = size(P, 2);
    correlations = zeros(p, 1);
    
    % Extracting upper triangle indices from a 160x160 matrix
    [i, j] = find(triu(ones(size(N)), 1));
    A1 = N(sub2ind(size(N), i, j));
    % Loop through columns of P
    for k = 1:p
        % Construct the Q_p distance matrix
        Q_p = abs(P(:,k) - P(:,k)');
        
        % Extract upper triangle elements from Q_p and correlate with N
        
        A2 = Q_p(sub2ind(size(Q_p), i, j));
        corrValue = corrcoef(A1, A2);
        correlations(k) = corrValue(2);
    end
    
    % Get the maximum correlation value
    maxCorr = max(correlations);
% end
