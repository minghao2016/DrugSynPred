%     clear
    addpath(genpath('code'));
    warning('off','all');    

    annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
    annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
    annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);

    [~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
    annotations.cellLines = annotations.cellLines(CL_perm, :);

        %% Read Monotherapy data
    %     fname = 'input/Dream/synergy/ch2_leaderBoard_monoTherapy.csv';
        fname = 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv';
    %     [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' ); %TODO: How to optimally combine replicates?

        [ Mono ] = read_singleMonoTherapy( annotations, fname );
        [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );
%%
%     [ interactome ] = readNetwork();
%     D2D = Construct_D2D(annotations, interactome);
%     [C2C, NodeWeights] = Construct_C2C(annotations, interactome, 'expression_only', true);    
%     Mono.EMax = fillinBlanks(Mono.EMax, C2C, D2D);
        
%     %%
%         synergy_threshold = 30;
%         X = Pair_synergy;
%         X(isinf(X)) = 0;
%         [Syn_pair_id, Syn_CL_id, Syn_vv] = find(X);
%         Syn_drug1 = Pairs(Syn_pair_id, 1);
%         Syn_drug2 = Pairs(Syn_pair_id, 2);
%         ck = Syn_CL_id;
%         Syn_labels = synergy_threshold <= Syn_vv;   
    %%
    Drug_Effect = cell2mat(Mono.EMax);
    Drug_Effect(isnan(Drug_Effect)) = 100;
    Drug_Effect = 100 - Drug_Effect; % effectiveness of drug in cell-lines
    
    [ Expr_mat, gene_names, expression_threshold ] = ReadExpr( 'Dream', annotations);
    
    Vert_compSL_Syn = ones(size(Pairs, 1), 1);
    for i = 1:size(Pairs, 1)
        [~, t1_celllines, t1_Effect] = find(Drug_Effect(Pairs(i, 1), :));        
        t1_targets = annotations.drugs.Target{Pairs(i, 1)};
        
        [t1_expr_mask, t1_expr_idx] = ismember(t1_targets, gene_names);     
        t1_expr_idx(~t1_expr_mask) = [];
        t1_targets(~t1_expr_mask) = [];
        
        
        [~, t2_celllines, t2_Effect] = find(Drug_Effect(Pairs(i, 2), :));        
        t2_targets = annotations.drugs.Target{Pairs(i, 2)};
        
        [t2_expr_mask, t2_expr_idx] = ismember(t2_targets, gene_names);         
        t2_expr_idx(~t2_expr_mask) = [];
        t2_targets(~t2_expr_mask) = [];
        
        t2_expr = Expr_mat(t2_expr_idx, t1_celllines);
        t1_expr = Expr_mat(t1_expr_idx, t2_celllines);



        idx1 = 0;
        p1 = ones(numel(t2_targets), 1);
        for k = 1:size(t2_expr, 1)
%             if(nnz(SL_net(SL_matches2(k), SL_matches1)) == 0)
%                 continue;
%             end
            idx1 = idx1 + 1;
            active = t2_expr(k, :) > expression_threshold;
            if(nnz(active) == 0 || nnz(~active) == 0)
                continue;
            end
                
            p1(idx1) = ranksum(t1_Effect(active), t1_Effect(~active), 'tail', 'left');
        end
        p1_sum = Edgington(p1);
        
        idx2 = 0;
        p2 = ones(numel(t1_targets), 1);
        for k = 1:size(t1_expr, 1)
%             if(nnz(SL_net(SL_matches1(k), SL_matches2)) == 0)
%                 continue;
%             end
            idx2 = idx2 + 1;
            active = t1_expr(k, :) > expression_threshold;
            if(nnz(active) == 0 || nnz(~active) == 0)
                continue;
            end
                
%             p2(idx2) = ranksum(t2_Effect(active), t2_Effect(~active), 'tail', 'left');
            p2(idx2) = ranksum(t2_Effect(active), t2_Effect(~active), 'tail', 'left');

        end
        p2_sum = Edgington(p2);
        
        pval = min(p1_sum, p2_sum);
%         pval2 = Edgington(x);
%         y = [sum(x), numel(x), pval, pval2]
%         pval = Edgington([p1; p2]);
                
        Vert_compSL_Syn(i, 1) = pval;        
    end    
    Vert_compSL_Syn(Vert_compSL_Syn < 0) = 1;
    
%%
    % Only works AFTER computing Synergy scores. If we move it to the
    % beginning results are horrible
    [ interactome ] = readNetwork();
    D2D = Construct_D2D(annotations, interactome);
    [C2C, NodeWeights] = Construct_C2C(annotations, interactome, 'expression_only', true);    
    X = cell2mat(fillinBlanks(Mono.EMax, C2C, D2D));
    
    X(isempty(X) | isnan(X) | isinf(X)) = 0; 
    X = normalize(X, 'dim', 1);    

    TT = Vert_compSL_Syn;
    TT(TT > 0.1) = 0;
    V{1} = full(spfun(@(x) -log10(x), TT));

%     V{1} = -log10(Vert_compSL_Syn);
    
    Z = zeros(size(Pairs, 1), size(annotations.cellLines, 1));
    for i = 1:size(Pairs, 1)
        for j = 1:size(annotations.cellLines, 1)
            Z(i, j) = V{1}(i) * max([X(Pairs(i, 1), j), X(Pairs(i, 2), j)], [], 2);
        end
    end    
%     Z = Modified_zscore(Z')';
%     Z(isnan(Z)) = 0;
%     Z = Z ./ max(nonzeros(Z));
    Z = 30*Z ./ prctile(nonzeros(Z), 95);   
    outPath = fullfile('output', 'predictions', 'leaderBoard');        
    exportResults( annotations, Pairs, Pair_names, Z, outPath, 'synergy_threshold', 20);
   