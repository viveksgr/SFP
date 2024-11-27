function M2 = SFP_restricter_hrf(M_main,M_ind,win)
if nargin<2
    win = 4:6;
end

[N, ~, V] = size(M_ind);
M2 = zeros(N, V);
fprintf('\n')
fprintf('Running HRF restriction...\n')


for n = 1:N
    % Odor responses for a given voxel
    restricted_values = squeeze(M_ind(n, win, :));
    [~, idx] = max(mean(restricted_values,2));
    actual_idx = idx + win(1)-1;
    M2(n, :) = squeeze(M_main(n, actual_idx, :));
end


fprintf('HRF restriction complete.\n')