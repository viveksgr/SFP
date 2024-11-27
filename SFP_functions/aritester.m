function [observed_ari,p_value,thresh]=aritester(idx,idx2);
observed_ari = rand_index(idx, idx2, 'adjusted');
num_permutations = 1000;
permuted_ari = zeros(num_permutations, 1);

for i = 1:num_permutations
    permuted_v1 = idx(randperm(length(idx)));
    permuted_ari(i) = rand_index(permuted_v1, idx2, 'adjusted');
end

thresh = prctile(permuted_ari,90);

p_value = sum(permuted_ari >= observed_ari) / num_permutations;