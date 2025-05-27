function [t_adj, p_adj, p_thresh] = fdr_t_scores(t_mat, dfs, q)
% FDR_T_SCORES  Convert row-wise t-scores to p-values, apply FDR correction,
%               and map the adjusted p-values back to t-scores.
%
%   Inputs:
%     t_mat  : [3 x T] matrix of t-scores
%     dfs    : [3 x 1] vector of degrees of freedom for each row
%     q      : desired FDR level (e.g. 0.05; default if omitted = 0.05)
%
%   Outputs:
%     t_adj   : [3 x T] matrix of t-scores corresponding to FDR-adjusted p-values
%     p_adj   : [3 x T] matrix of FDR-adjusted p-values
%     p_thresh: scalar p-value threshold returned by fdr_benjhoc
%
%   Requires the function fdr_benjhoc (as provided).

    if nargin<3 || isempty(q)
        q = 0.05;
    end

    [nRows, nT] = size(t_mat);
    if nRows~=numel(dfs)
        error('Length of dfs must match number of rows in t_mat');
    end

    % 1) Convert t-scores to two-tailed p-values
    p_mat = nan(size(t_mat));
    for r = 1:nRows
        df = dfs(r);
        t_vec = t_mat(r,:);
        p_mat(r,:) = 2 * (1 - tcdf(abs(t_vec), df));
    end

    % 2) Flatten p-values and FDR-correct
    p_vec = p_mat(:);
    [p_thresh, p_correct] = fdr_benjhoc(p_vec, q);

    % 3) Reshape adjusted p-values back to matrix
    p_adj = reshape(p_correct, [nRows, nT]);

    % 4) Convert adjusted p-values back to t-scores,
    %    preserving sign of original t_mat
    t_adj = nan(size(t_mat));
    for r = 1:nRows
        df = dfs(r);
        sign_vec = sign(t_mat(r,:));
        % two-tailed: p = 2*(1 - tcdf(|t|,df))  =>  |t| = tinv(1 - p/2, df)
        t_mag = tinv(1 - p_adj(r,:)/2, df);
        t_adj(r,:) = sign_vec .* t_mag;
    end
end
