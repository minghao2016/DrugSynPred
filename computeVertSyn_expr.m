function [ VertSyn_Expr, VertSyn_Expr_Rest ] = computeVertSyn_expr( annotations, Pairs )
    [ ACSN ] = import_ACSN();
    [ transcriptional_gene_signature, LINCS_genes ] = computeTransSig( annotations );
    ACSN_genes = [ACSN.class_genes{:}];
    [mask, idx]  = ismember(ACSN_genes, LINCS_genes);
    idx(~mask) = [];
    ACSN_genes(~mask) = [];

    %%
    respExpr_mat = zeros(numel(LINCS_genes), size(annotations.drugs, 1));
    respExprRest_mat = zeros(numel(ACSN_genes), size(annotations.drugs, 1));
    for i = 1:size(annotations.drugs, 1)
        if(isempty(transcriptional_gene_signature{i, 24}))
            continue;
        end
        respExpr_mat(:, i) = transcriptional_gene_signature{i, 24};
        respExprRest_mat(:, i) = transcriptional_gene_signature{i, 24}(idx);
    end
    %% For all genes    
    [CC, pval] = corr(respExpr_mat);
    CC( abs(CC - 1) < eps | isnan(CC) | pval > 1e-3 | CC <= 0) = 0;
        
    VertSyn_Expr = zeros(size(Pairs, 1), 1);
    for i = 1:size(Pairs, 1)        
        VertSyn_Expr(i) = CC(Pairs(i, 1), Pairs(i, 2));        
    end

    %% For the set of genes restricted to ACSN classes
    [CC, pval] = corr(respExprRest_mat);
    CC( isnan(CC) ) = 0;
    CC( abs(CC - 1) < eps | pval > 1e-3 | CC <= 0) = 0;
        
    VertSyn_Expr_Rest = zeros(size(Pairs, 1), 1);
    for i = 1:size(Pairs, 1)        
        VertSyn_Expr_Rest(i) = CC(Pairs(i, 1), Pairs(i, 2));        
    end    
    
%     [CC, pval] = corr(VertSyn_Expr_Rest, [sum(Pair_synergy > 20, 2) ./ sum(~isinf(Pair_synergy), 2), sum(Pair_synergy > 30, 2) ./ sum(~isinf(Pair_synergy), 2) , sum(Pair_synergy > 40, 2) ./ sum(~isinf(Pair_synergy), 2)])
end


