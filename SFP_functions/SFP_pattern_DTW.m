function DTWout   = SFP_pattern_DTW(FMat, G,optout)
if nargin<3
    optout = 1000;
end
% Ensure G is a column vector
if size(G, 2) > 1
    G = G';
end

% Identify the unique conditions
conditions = unique(G);
NCond = length(conditions);

% Initialize output
DTWout = cell(NCond, 1);
% actual_statistic = zeros(NCond, 1);
% permuted_stats = zeros(NCond, optout);

[~, T] = size(FMat);  % Extracting the number of time points

for i = 1:NCond
    current_condition = conditions(i);

    % T1: Subset of FMat for the current condition
    T1 = FMat(G == current_condition, :);

    % T2: Average time traces for all conditions except the current one
    T2 = zeros(NCond-1, T);  % Pre-allocating space
    idx = 1;  % Index to keep track of the current row for T2

    for j = 1:NCond
        if j ~= i
            avg_trace = mean(FMat(G == conditions(j), :), 1);
            T2(idx, :) = avg_trace;
            idx = idx + 1;
        end
    end

    % Compute DTW
    % DTWout(i) = SFP_computeDTW_bootstrap(T1, T2,optout);
    % [DTWout(i), actual_statistic(i), permuted_stats(i,:)] = SFP_computeDTW_perm(T1, T2, optout, 1000, false);
    DTWout{i} = SFP_computeDTW(T1, T2, optout,false);
end
end
