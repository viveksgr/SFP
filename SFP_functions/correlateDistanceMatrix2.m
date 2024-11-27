function maxCorr = correlateDistanceMatrix2(N, P)
    % N: 4320x4320 correlation matrix
    % P: 4320xp matrix

    p = size(P, 2);
    correlations = zeros(p, 1);
    
    % Extracting upper triangle indices from a 4320x4320 matrix
    [i, j] = find(triu(ones(size(N)), 1));
    linearIndices = sub2ind(size(N), i, j);
    nUpperTri = N(linearIndices);
    
    % Loop through columns of P
    for k = 1:p
        % Construct the Q_p distance matrix and extract upper triangle
        qUpperTri = abs(P(:,k) - P(:,k)');
        qUpperTri = qUpperTri(linearIndices);
        
        % Compute correlation
        correlations(k) = corr(nUpperTri, qUpperTri);
    end
    
    % Get the maximum correlation value
    maxCorr = max(correlations);
end
