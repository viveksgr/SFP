function cooccurrence = sfp_cooccurrenceMatrix(idx, idx2)
    labels1 = unique(idx);
    labels2 = unique(idx2);
    
    cooccurrence = zeros(numel(labels1), numel(labels2));
    
    for i = 1:numel(labels1)
        for j = 1:numel(labels2)
            cooccurrence(i, j) = sum((idx == labels1(i)) & (idx2 == labels2(j)));
        end
    end
end