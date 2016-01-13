function [ VertSyn ] = computeVertSyn_topo( annotations, interactome, Pairs )
    % combine both targets and their functional neighborhood
    drug_target_net_indices = arrayfun(@(x) find(ismember(interactome.vertex_names, annotations.drugs.Target{x}))', 1:size(annotations.drugs, 1), 'UniformOutput', false);
    Expanded_targets = arrayfun(@(d) unique([[interactome.Pruned_neighbors{drug_target_net_indices{d}}], drug_target_net_indices{d}]), 1:size(annotations.drugs, 1), 'UniformOutput', false)'; 

    % compute vertical synergy scores (if drugs target upstream/downstream
    % elements in network)
    VertSyn = zeros(size(Pairs, 1), 1);
    for i = 1:size(Pairs, 1)
        d1 = Pairs(i, 1);
        d2 = Pairs(i, 2);
        common = numel(intersect(Expanded_targets{d1}, Expanded_targets{d2}));
        if(common == 0)
            continue;
        end
        logPval = -log10(hygecdf(common-1, numel(interactome.vertex_names), numel(Expanded_targets{d1}), numel(Expanded_targets{d2}), 'upper'));
        VertSyn(i) = logPval;
    end    
end


