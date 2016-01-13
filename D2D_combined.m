function [ D2D ] = D2D_combined( annotations, interactome )
    stitch = load('input/STITCH/drug_comb_mapped');
    stitchD2D = sparse(stitch(:,1), stitch(:,2), stitch(:, 3), 119, 119);
    stitchD2D = max(stitchD2D, stitchD2D');
    stitchD2D = full(stitchD2D ./ max(stitchD2D(:)));
    stitchD2D(~stitchD2D) = nan;
    stitchD2D = stitchD2D ./ nanmax(nonzeros(stitchD2D));

    %%
%     drug_target_net_indices = arrayfun(@(x) find(ismember(interactome.vertex_names, annotations.drugs.Target{x}))', 1:size(annotations.drugs, 1), 'UniformOutput', false);
%     Expanded_targets = arrayfun(@(d) unique([[interactome.Pruned_neighbors{drug_target_net_indices{d}}], drug_target_net_indices{d}]), 1:size(annotations.drugs, 1), 'UniformOutput', false)'; 
%     
%     targetD2D = nan(size(annotations.drugs, 1));
%     for i = 1:size(annotations.drugs, 1)
%         for j = i+1:size(annotations.drugs, 1)
%             common = numel(intersect(Expanded_targets{i}, Expanded_targets{j}));
%             logPval = -log10(hygecdf(common-1, numel(interactome.vertex_names), numel(Expanded_targets{i}), numel(Expanded_targets{j}), 'upper'));
%             targetD2D(i, j) = logPval;
%         end
%     end
%     targetD2D = max(targetD2D, targetD2D');
%     targetD2D = targetD2D ./ nanmax(nonzeros(targetD2D));
    
%%
    % Compute RWR over interactome from targets to identify topological signature of each drug
        par.alpha = 0.9;
        n = size(interactome.A, 1);
        P = interactome.A'*spdiags(spfun(@(x) 1./x, sum(interactome.A, 2)), 0, n, n);

        % Handling the dangling nodes ...
        e_T = ones(1, n);    
        d_T = e_T - e_T*P;
        P = P + diag(d_T);
        Q = (1-par.alpha).*inv(eye(n) - par.alpha*P);

        topological_signatures = zeros(n, size(annotations.drugs, 1));
        for i = 1:size(annotations.drugs, 1)
            primary_targets = annotations.drugs.Target{i};
            [~, src_nodes] = ismember(primary_targets, interactome.vertex_genes);
            src_nodes(~src_nodes) = [];
            if(numel(src_nodes) == 0)
                continue;
            end
            e_src = sparse(src_nodes, 1, 1, n, 1); e_src = e_src ./ sum(e_src);
            topological_signatures(:, i) = Q*e_src;
        end
        TT = log10(topological_signatures);
        TT(isinf(TT)) = 0;
        TT(TT ~= 0) = TT(TT ~= 0) - min(nonzeros(TT));
        targetD2D = partialcorr(TT, mean(TT, 2));
        
%         targetD2D = partialcorr(topological_signatures, mean(topological_signatures, 2));
        pos_targetD2D = targetD2D;
        pos_targetD2D(pos_targetD2D < 0) = nan;
        pos_targetD2D = pos_targetD2D - diag(diag(pos_targetD2D));
        pos_targetD2D = pos_targetD2D ./ nanmax(nonzeros(pos_targetD2D));
    %%
    classD2D = nan(size(annotations.drugs, 1));
    drug_classes = cellfun(@(c) strsplit(c, ';'), annotations.drugs.MoA, 'UniformOutput', false);
    classes = {'ALK','ANTIMETABOLITE','AUTO','HDAC','PROT','TOPO','TUBB'};

    
    for i = 1:numel(classes)
        rows = find(cellfun(@(c) ~isempty(find(ismember(c, classes{i}))), drug_classes));
        classD2D(rows, rows) = 1;        
    end    
    
    %%
    D2D = reshape(nanmean([stitchD2D(:), pos_targetD2D(:), classD2D(:)], 2), size(stitchD2D, 1), size(stitchD2D, 1));    
    D2D(isnan(D2D)) = 0;

end

