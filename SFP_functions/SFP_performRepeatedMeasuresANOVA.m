function pValue = SFP_performRepeatedMeasuresANOVA(cellArray)
    % Validate input
    if ~iscell(cellArray) || size(cellArray, 1) ~= 3 || size(cellArray, 2) ~= 5
        error('Input must be a 3x5 cell array.');
    end
    
    % Prepare the data
    [nSubjects, nConditions] = size(cellArray);
    subjectID = [];
    conditionID = [];
    values = [];

    for i = 1:nSubjects
        for j = 1:nConditions
            currentLength = length(cellArray{i, j});
            subjectID = [subjectID; repmat(i, currentLength, 1)];
            conditionID = [conditionID; repmat(j, currentLength, 1)];
            values = [values; cellArray{i, j}];
        end
    end

    % Create a table for the mixed-effects model
    tbl = table(subjectID, conditionID, values, 'VariableNames', {'Subject', 'Condition', 'Values'});
    tbl.Subject = nominal(tbl.Subject);  % Treat subject ID as a categorical variable
    tbl.Condition = nominal(tbl.Condition);  % Treat condition ID as a categorical variable

    % Fit the mixed-effects model
    lme = fitlme(tbl, 'Values ~ Condition + (1|Subject)');

    % Perform ANOVA to test for the effect of Condition
    results = anova(lme);
    pValue = results.pValue(2);  % p-value for the condition effect

    FStatistic = results.FStat(2);

    % Access degrees of freedom associated with the F-statistic
    DF1 = results.DF1(2); % Degrees of freedom for the numerator (between groups)
    DF2 = results.DF2(2); % Degrees of freedom for the denominator (within groups, error)

    % Display or format the result for reporting
    fprintf('The F-statistic for the condition effect is F(%d, %d) = %.3f, p = %.3f.\n', ...
        DF1, DF2, FStatistic, pValue);

    posthocResults = performPosthocComparisonsFDR(lme);
    % Display results
    fprintf('ANOVA p-value for Condition effect: %f\n', pValue);

    % Generate box plots
    figure;
    boxplot(tbl.Values, {tbl.Condition, tbl.Subject}, 'factorgap', [5, 2], 'factorseparator', [1]);
    title('Box Plots of Values by Condition and Subject');
    xlabel('Condition and Subject');
    ylabel('Values');
    grid on;
    
    return
end


function posthocResults = performPosthocComparisons(lme)
    % Extract the fixed effects estimates and their covariance matrix
    [beta,~,stats] = fixedEffects(lme);
    
    % Extract standard errors and the degrees of freedom
    SE = stats.SE;
    df = stats.DF;

    % Number of conditions
    numConditions = length(beta);

    % Pre-allocate matrix to hold p-values and confidence intervals
    comparisons = nchoosek(1:numConditions, 2); % All pairwise comparisons
    numComparisons = size(comparisons, 1);
    posthocResults = table(comparisons, zeros(numComparisons, 1), zeros(numComparisons, 2), ...
        'VariableNames', {'Comparison', 'pValue', 'ConfidenceInterval'});

    % Perform each pairwise comparison
    for i = 1:numComparisons
        c1 = comparisons(i, 1);
        c2 = comparisons(i, 2);
        
        % Difference in estimates
        diffEst = beta(c1) - beta(c2);
        
        % Standard error of the difference
        SE_diff = sqrt(SE(c1)^2 + SE(c2)^2);
        
        % t statistic
        tStat = diffEst / SE_diff;
        
        % p-value (two-tailed)
        pValue = 2 * tcdf(-abs(tStat), df(1));
        adjustedP = min(pValue * numComparisons, 1);
        
        % 95% Confidence Interval
        CI = [diffEst - 1.96 * SE_diff, diffEst + 1.96 * SE_diff];
        
        % Store results
        posthocResults.pValue(i) = adjustedP;
        posthocResults.ConfidenceInterval(i, :) = CI;
    end

    % Display results
    disp(posthocResults);
end


function posthocResults = performPosthocComparisonsFDR(lme)
    % Extract fixed effects estimates, covariance matrix, standard errors, and degrees of freedom
    [beta,n,stats] = fixedEffects(lme);
    SE = stats.SE;
    df = stats.DF(1);  % Assuming consistent DF for simplicity

    % All pairwise comparisons
    comparisons = nchoosek(1:length(beta), 2);
    pValues = zeros(size(comparisons, 1), 1);

    % Calculate p-values for each comparison
    for i = 1:size(comparisons, 1)
        diffEst = beta(comparisons(i, 1)) - beta(comparisons(i, 2));
        SE_diff = sqrt(SE(comparisons(i, 1))^2 + SE(comparisons(i, 2))^2);
        tStat = diffEst / SE_diff;
        pValues(i) = 2 * tcdf(-abs(tStat), df);
    end

    % % Adjust p-values for FDR
     % adjustedPValues = pValues;
    adjustedPValues = mafdr(pValues, 'BHFDR', true);

    % Store results
    posthocResults = table(comparisons, adjustedPValues, 'VariableNames', {'Comparisons', 'AdjustedPValues'});
    disp(posthocResults);
end
