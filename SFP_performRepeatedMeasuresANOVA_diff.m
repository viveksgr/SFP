function [adjustedPValues, t_st,df] = SFP_performRepeatedMeasuresANOVA_diff(cellArray)

% Example to convert cell array to table
numSubjects = size(cellArray,1);
numAreas = size(cellArray,2);
numConditions = size(cellArray,3); % Q1, Q2, Q3
measurements = [];
subjectID = [];
areaID = [];
conditionID = [];

for s = 1:numSubjects
    for a = 1:numAreas
        for c = 1:numConditions
            currentData = cellArray{s, a, c};
            numDataPoints = length(currentData);
            
            measurements = [measurements; currentData];
            subjectID = [subjectID; repmat(s, numDataPoints, 1)];
            areaID = [areaID; repmat(a, numDataPoints, 1)];
            conditionID = [conditionID; repmat(c, numDataPoints, 1)];
        end
    end
end

tbl = table(subjectID, areaID, conditionID, measurements, ...
            'VariableNames', {'Subject', 'AnatomicalArea', 'Condition', 'Measurement'});
tbl.Subject = nominal(tbl.Subject);
tbl.AnatomicalArea = nominal(tbl.AnatomicalArea);
tbl.Condition = nominal(tbl.Condition);

% If Condition is stored as numeric identifiers and you want them categorical
tbl.Condition = categorical(tbl.Condition, {'1', '2', '3'}, 'Ordinal', false);
tbl.Condition = nominal(tbl.Condition);
% Set '3' as the reference category
tbl.Condition = reorderlevels(tbl.Condition, {'3', '1', '2'});


areas = unique(tbl.AnatomicalArea);  % Get unique anatomical areas
pValues = zeros(length(areas), 2);  % Store p-values from each ANOVA
t_st = zeros(length(areas), 2); 
df = zeros(length(areas), 1); 
for i = 1:length(areas)
    areaData = tbl(tbl.AnatomicalArea == areas(i), :);  % Data for current area
    
    % Fit a linear mixed-effects model (or use fitrm for repeated measures)
    lme = fitlme(areaData, 'Measurement ~ Condition + (1|Subject)');
    
    pValues(i,1)=lme.Coefficients{2,6};
    pValues(i,2)=lme.Coefficients{3,6};
    t_st(i,1)=lme.Coefficients{2,4};
    t_st(i,2)=lme.Coefficients{3,4};
    df(i)=lme.Coefficients{2,5};
end

% Apply Bonferroni correction
adjustedPValues = min(pValues *2* length(areas), 1);

end
