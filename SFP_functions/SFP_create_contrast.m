function contrastMatrix = SFP_create_contrast(Q,numAreas)
% Contrast matrix
contrastMatrix = zeros(1, numAreas * 3); % Initialize for 3 conditions across 5 areas

% Set the contrasts
for i = 1:numAreas
    indexQ1 = 1 + (i - 1) * 3; % Adjust index for each block of condition terms
    indexQ2 = 2 + (i - 1) * 3;
    indexQ3 = 3 + (i - 1) * 3;
    
    contrastMatrix(indexQ1) = Q(1); % Coefficient for Q1 in each area
    contrastMatrix(indexQ2) = Q(2);  % Coefficient for Q2 in each area
    contrastMatrix(indexQ3) = Q(3); % Coefficient for Q3 in each area
end
end
