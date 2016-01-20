function [ Horiz_compSL_Syn ] = computeHorSyn_DrugEffect( annotations, Pairs, Drug_Effect )
    [ Expr_mat, gene_names, expression_threshold ] = ReadExpr( 'Dream', annotations);    

    empty_cols = find(sum(Expr_mat) == 0);
    Horiz_compSL_Syn = ones(size(Pairs, 1), 1);
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



        if(0 < numel(t1_targets))
            idx1 = 0;
            p1 = ones(numel(t2_targets), 1);
            for k = 1:size(t2_expr, 1)
    %             if(nnz(SL_net(SL_matches2(k), SL_matches1)) == 0)
    %                 continue;
    %             end
                min_pval = 1;
                [~, perm] = sort(t2_expr(k, :), 'descend');
                cuts = nnz(t2_expr(k, :)> expression_threshold); 
                if(cuts == 0 || cuts == 85)
                    continue;
                end


%                 for cuts = 1:numel(perm)-empty_cols-1
                    Eff_active = t1_Effect(perm(1:cuts));
                    Eff_inactive = t1_Effect(perm(cuts+1:end));
                    curr_pval = ranksum(Eff_active, Eff_inactive, 'tail', 'left');
                    if(curr_pval < min_pval)
                        min_pval = curr_pval;                    
                    end
%                 end
                idx1 = idx1 + 1;               
                p1(idx1) = min_pval;
            end
            p1_sum = Edgington(p1);
        else
            p1_sum = 1;
        end

        if(0 < numel(t2_targets))
            idx2 = 0;
            p2 = ones(numel(t1_targets), 1);
            for k = 1:size(t1_expr, 1)
                min_pval = 1;
                [~, perm] = sort(t1_expr(k, :), 'descend');
                cuts = nnz(t1_expr(k, :) > expression_threshold);
                if(cuts == 0 || cuts == 85)
                    continue;
                end
%                 for cuts = 1:numel(perm)-empty_cols-1
                    Eff_active = t2_Effect(perm(1:cuts));
                    Eff_inactive = t2_Effect(perm(cuts+1:end));
                    curr_pval = ranksum(Eff_active, Eff_inactive, 'tail', 'left');
                    if(curr_pval < min_pval)
                        min_pval = curr_pval;                    
                    end
%                 end
                idx2 = idx2 + 1;               
                p2(idx2) = min_pval;
            end
            p2_sum = Edgington(p2);
        else
            p2_sum = 1;
        end
        
        pval = min(p1_sum, p2_sum);
%         pval2 = Edgington(x);
%         y = [sum(x), numel(x), pval, pval2]
%         pval = Edgington([p1; p2]);
                
        Horiz_compSL_Syn(i, 1) = pval;        
    end    
    Horiz_compSL_Syn(Horiz_compSL_Syn < 0) = 1;

end

