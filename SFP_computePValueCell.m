function pValue = SFP_computePValue(M1, M2, df)
    N = length(M1);  % Number of regression results
    betaMean = mean(M1);  % Mean beta weight
    SE = M1 ./ M2;  % Calculate individual standard errors
    SE_MeanBeta = sqrt(sum(SE.^2) / N) / sqrt(N);  % SE of the mean beta weight
    tScore = betaMean / SE_MeanBeta;  % T-score for the mean beta weight
    pValue = 2 * (1 - tcdf(abs(tScore), df));  % P-value from t-distribution
end