function crosspval = SFP_crosspvalcalc(weights,tscores,df)

% Extract weights and standard errors for x1 and x2
weight_x1 = weights(1);
weight_x2 = weights(2);

se_x1 = weights(1)/tscores(1);
se_x2 =  weights(2)/tscores(2);

% Calculate the difference between the weights
weight_diff = weight_x1 - weight_x2;

% Calculate the standard error of the difference
se_diff = sqrt(se_x1^2 + se_x2^2);

% Calculate the t-score for the difference
t_score_diff = weight_diff / se_diff;

% Calculate the p-value for the two-tailed t-test
crosspval = 2 * (1 - tcdf(abs(t_score_diff), df));
