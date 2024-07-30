function [testStatistic, combinedPValue] = fishersCombinedTest(pvalues)
    % Check that the input is a 36x1 vector
    % Calculate Fisher's combined test statistic
    
    testStatistic = -2 * sum(log(pvalues));

    % Degrees of freedom is twice the number of p-values
    df = 2 * numel(pvalues);

    % Calculate the combined p-value using the chi-square distribution
    combinedPValue = 1 - chi2cdf(testStatistic, df);
end
