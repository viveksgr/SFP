function p_value = ARC_computePValueOneTailed(r_eff, n, num_samples)
    % Chance level accuracy
    chance_accuracy = 1/n;
    
    % Assuming binomial distribution under the null hypothesis
    % Variance of the binomial distribution is p*(1-p)/N, where p is the chance accuracy,
    % and N is the total number of trials (from num_samples)
    variance = chance_accuracy * (1 - chance_accuracy) / num_samples;
    
    % Standard deviation of accuracies
    std_dev = sqrt(variance);
    
    % Compute the z-score
    z_score = (r_eff - chance_accuracy) / std_dev;
    
    % Convert the z-score to a p-value for a one-tailed test
    p_value = 2*(1 - normcdf(abs(z_score), 0, 1)); % Only interested in performance above chance
end
