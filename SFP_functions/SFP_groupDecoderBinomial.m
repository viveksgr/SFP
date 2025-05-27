function [Zgroup,pGroup,zSubs] = SFP_groupDecoderBinomial(acc)
% acc      : vector of per-subject accuracies (decoding proportion)
% nTrials  : vector of number of trials per subject
% nCond    : vector of number of conditions (chance = 1/nCond) per subject
% zSubs    : subject-level z-scores (one-tailed, H0 : accuracy = chance)
% Zgroup   : Stouffer weighted Z (weights = sqrt(nTrials))
% pGroup   : group one-tailed p-value
% acc = [0.0195 0.0238 0.0435 0.125 0.135 0.155];
    nTrials  = [4560,4320,4320,     200,200,200,200];
    nCond    = [160,160,160,        10,10,10,10];

    pChance   = 1 ./ nCond;                 % chance level for each subject
    varChance = pChance .* (1-pChance) ./ nTrials;
    zSubs     = (acc - pChance) ./ sqrt(varChance);   % subject z

    w         = sqrt(nTrials);              % Stouffer weights
    Zgroup    = sum(w .* zSubs) / sqrt(sum(w.^2));
    pGroup    = 1 - normcdf(Zgroup);        % upper-tail p
end
