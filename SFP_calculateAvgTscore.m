function avgTscore = SFP_calculateAvgTscore(M1, M2)
    % Calculate the mean beta weight
    betaMean = mean(M1(:));
    
    % Calculate individual standard errors from beta and t-score matrices
    SE_individual = M1 ./ M2;
    
    % Calculate the RMS of the standard errors
    SE_RMS = sqrt(median(SE_individual(:).^2,"omitnan"));
    
    % Calculate the t-score for the mean beta weight using the RMS of SEs
    avgTscore = betaMean / SE_RMS;
    
    % Display the calculated average t-score
    % fprintf('The average t-score for the mean beta weight using RMS of SEs is: %f\n', avgTscore);
end
