function [behavnew, behavAvg] = averageBehavWithUnionOdors(behav)
% AVERAGEBEHAVWITHUNIONODORS  Averages behavioral ratings across multiple subjects
% by taking the union of all odors and descriptors used by any subject.
%
%   behavAvg = AVERAGEBEHAVWITHUNIONODORS(behav)
%
%   INPUT:
%       behav : A structure array, one element per subject, each containing:
%               .cid        -> [N x 1] array of odor IDs for that subject
%               .percept    -> [M x 1] array of descriptor IDs or names
%               .ratings    -> [N x M] rating matrix (odor vs descriptor)
%
%   OUTPUT:
%       behavAvg : A structure with fields:
%                   .cid        -> Union of all odor IDs across subjects
%                   .percept    -> Union of all descriptors across subjects
%                   .ratings    -> [numOdors x numDescripts] matrix of average ratings
%                                  (NaN if an odor–descriptor pair is absent for all subjects).
%
%   EXAMPLE:
%       % Suppose behav is a 1x3 structure for 3 subjects, each with .cid, .percept,
%       % and .ratings, but they may have partially overlapping sets of odors and descriptors.
%       behavAvg = averageBehavWithUnionOdors(behav);
%       disp(behavAvg);
%

    numSubjects = numel(behav);
    if numSubjects < 1
        error('The behav structure must contain at least one subject.');
    end

    % 1) Build the union of all odor IDs
    allOdors = [];
    for s = 1:numSubjects
        allOdors = [allOdors; behav(s).cid];
    end
    % Convert to column, find unique, preserve order of first occurrence
    allOdors = unique(allOdors, 'stable');
    behavAvg.cid = allOdors;  % union of odors in stable order
    
    % 2) Build the union of all descriptors
    allDescriptors = [];
    for s = 1:numSubjects
        allDescriptors = [allDescriptors; behav(s).percepts];
    end
    allDescriptors = unique(allDescriptors, 'stable');
    behavAvg.percept = allDescriptors;

    numOdors = length(allOdors);
    numDescripts = length(allDescriptors);

    % 3) Initialize the average ratings matrix with NaNs
    behavAvg.ratings = NaN(numOdors, numDescripts);

    % 4) Accumulate ratings by looping over each subject
    %    We'll use sums and counts to do averaging only after collecting from all subjects
    sumRatings = zeros(numOdors, numDescripts);
    countRatings = zeros(numOdors, numDescripts);

    for s = 1:numSubjects
        cidsSubj = behav(s).cid;        % e.g. Nx1
        percSubj = behav(s).percepts;    % e.g. Mx1
        ratingsSubj = behav(s).ratings; % [N x M]

        % For each odor-descriptor pair in this subject, map to row/col in big matrix
        for iO = 1:length(cidsSubj)
            odorID = cidsSubj(iO);
            rowIdx = find(allOdors == odorID, 1, 'first');
            if isempty(rowIdx), continue; end

            for iP = 1:length(percSubj)
                descID = percSubj(iP);
                colIdx = find(strcmp(allDescriptors, descID), 1, 'first');
                if isempty(colIdx), continue; end

                val = ratingsSubj(iO, iP);
                % Update sum and count for averaging
                if ~isnan(val)
                    sumRatings(rowIdx, colIdx)   = sumRatings(rowIdx, colIdx) + val;
                    countRatings(rowIdx, colIdx) = countRatings(rowIdx, colIdx) + 1;
                end
            end
        end
    end

    % 5) Compute average where count is > 0
    nonzeroMask = (countRatings > 0);
    behavAvg.ratings(nonzeroMask) = sumRatings(nonzeroMask) ./ countRatings(nonzeroMask);
    
    behavnew = rebuildBehavFromAvg(behav, behavAvg);
end

function newBehav = rebuildBehavFromAvg(oldBehav, behavAvg)
% REBUILDBEHAVFROMAVG Rebuilds a new behavior structure array using averaged data
%                     while preserving each subject's original odor/descriptor order.
%
%   newBehav = rebuildBehavFromAvg(oldBehav, behavAvg)
%
%   INPUTS:
%       oldBehav : A 1xS structure array (original format), each entry having:
%                     .cid       -> [N x 1] array of odor IDs for that subject
%                     .percept   -> [M x 1] array of descriptor IDs or names
%                     .ratings   -> [N x M] rating matrix (unused here except for shape)
%       behavAvg : A structure with fields:
%                     .cid       -> Union of odor IDs (e.g. length = 195),
%                     .percept   -> Union of descriptor IDs,
%                     .ratings   -> [#odors x #descriptors] matrix of averaged ratings.
%
%   OUTPUT:
%       newBehav : A 1xS structure with the same shape as oldBehav, but each subject's
%                  .ratings now holds data extracted from behavAvg, preserving the
%                  subject's original row/column ordering in .cid and .percept.
%                  If an odor or descriptor isn't found in behavAvg, the rating is NaN.
%
%   EXAMPLE:
%       % Suppose oldBehav is the original data for 3 subjects
%       % and behavAvg is the union-based average structure.
%       newBehav = rebuildBehavFromAvg(oldBehav, behavAvg);
%       disp(newBehav);

    numSubjects = numel(oldBehav);
    newBehav(1:numSubjects) = struct('cid', [], 'percepts', [], 'ratings', []);

    for s = 1:numSubjects
        % Original odor IDs and descriptors for subject s
        oldCid = oldBehav(s).cid;         % e.g. 160x1
        oldPercept = oldBehav(s).percepts; % e.g. Nx1

        nOdors = length(oldCid);
        nDescs = length(oldPercept);

        % Initialize ratings with NaN
        newRatings = NaN(nOdors, nDescs);

        % Loop over each odor for subject s
        for i = 1:nOdors
            odorID = oldCid(i);

            % Locate this odor in behavAvg.cid (assumes numeric odor IDs)
            rowIdx = find(behavAvg.cid == odorID, 1);
            if isempty(rowIdx)
                % This odor is not present in behavAvg -> remain NaN
                continue;
            end

            % Loop over each descriptor for subject s
            for j = 1:nDescs
                descID = oldPercept(j);

                % Locate this descriptor in behavAvg.percept (assumes numeric)
                colIdx = find(strcmp(behavAvg.percepts,descID), 1, 'first');

                if isempty(colIdx)
                    % This descriptor is not present in behavAvg -> remain NaN
                    continue;
                end

                % Extract the averaged rating
                newRatings(i, j) = behavAvg.ratings(rowIdx, colIdx);
            end
        end

        % Build the new structure entry for subject s
        newBehav(s).cid = oldCid;               % preserve original odor order
        newBehav(s).percepts = oldPercept;       % preserve original descriptor order
        newBehav(s).ratings = newRatings;       % newly filled from behavAvg
    end
end

