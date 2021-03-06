    clear
    addpath(genpath('code'));
    warning('off','all');    

    annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
    annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
    annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);

    [~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
    annotations.cellLines = annotations.cellLines(CL_perm, :);


    load datasets

%%
    fname = 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv';
    [ Mono ] = read_singleMonoTherapy( annotations, fname );

    %%
    [ VertSyn ] = computeVertSyn_topo( annotations, interactome, Pairs );
    
    %
    
%% Enrichment of synthetic lethal interactions among functional/physical neighborhood of drug targets (parallel pathways)
    fprintf('Computing Synergy for drug pairs based on SL interactions\n');
%     fd = fopen('input/SynLethDB/sl_human', 'r');
%     C = textscan(fd, '%s %s');
%     fclose(fd);

    [Table, header] = my_dlmread('input/SynLethDB/sl_human');
    SL_src = Table(:, 1);
    SL_dst = Table(:, 3);
    SL_conf = cellfun(@(x) str2double(x), Table(:, 9));

    SL_conf_threshold = 0.2;
    SL_src(SL_conf < SL_conf_threshold) = [];
    SL_dst(SL_conf < SL_conf_threshold) = [];
    SL_conf(SL_conf < SL_conf_threshold) = [];
    
    
    SL_nodes = union(SL_src, SL_dst);
    [~, ii] = ismember(SL_src, SL_nodes);
    [~, jj] = ismember(SL_dst, SL_nodes);
    SL_net = sparse(ii, jj, 1, numel(SL_nodes), numel(SL_nodes));
    SL_net = max(SL_net, SL_net');
    
    population_size = nchoosek(numel(SL_nodes), 2);
    total_success = nnz(SL_net)/2;

    Vert_syn = zeros(size(Pairs, 1), 1);
    Horiz_syn = zeros(size(Pairs, 1), 1);    
    for i = 1:size(Pairs, 1)
        N1_names = interactome.vertex_names(Drug_target_neighbors{Pairs(i, 1)});
        [~, SL_matches1] = ismember(N1_names, SL_nodes);
        SL_matches1(~SL_matches1) = [];
        if(isempty(SL_matches1))
            continue;
        end
        
        N2_names = interactome.vertex_names(Drug_target_neighbors{Pairs(i, 2)});
        [~, SL_matches2] = ismember(N2_names, SL_nodes);
        SL_matches2(~SL_matches2) = [];
        if(isempty(SL_matches2))
            continue;
        end        
        
        success_size = numel(intersect(N1_names, N2_names));
        Horiz_syn(i) = -log10(hygecdf(success_size-1, numel(interactome.vertex_names), numel(N1_names), numel(N2_names), 'upper'));                

        success_size = nnz(SL_net(SL_matches1, SL_matches2));
        Vert_syn(i) = -log10(hygecdf(success_size-1, nchoosek(numel(SL_nodes), 2), nnz(SL_net)/2, numel(SL_matches1)*numel(SL_matches2), 'upper'));                
    end

    
%%
    
    Drug_Effect = cell2mat(Mono.EMax);
    Drug_Effect(isnan(Drug_Effect)) = 100;
    Drug_Effect = 100 - Drug_Effect; % effectiveness of drug in cell-lines


    expression_threshold = 5.5;
    
    Vert_Syn_Confirm_pval = ones(size(Pairs, 1), 1);
    for i = 1:size(Pairs, 1)
        [~, t1_celllines, t1_Effect] = find(Drug_Effect(Pairs(i, 1), :));
        
%         t1_targets = interactome.vertex_names(Drug_target_neighbors{Pairs(i, 1)});
        t1_targets = annotations.drugs.Target{Pairs(i, 1)};
        
        [t1_expr_mask, t1_expr_idx] = ismember(t1_targets, Dream_gene_names);     
        t1_expr_idx(~t1_expr_mask) = [];
        t1_targets(~t1_expr_mask) = [];
        
        
        [~, t2_celllines, t2_Effect] = find(Drug_Effect(Pairs(i, 2), :));
        
%         t2_targets = interactome.vertex_names(Drug_target_neighbors{Pairs(i, 2)});
        t2_targets = annotations.drugs.Target{Pairs(i, 2)};
        
        [t2_expr_mask, t2_expr_idx] = ismember(t2_targets, Dream_gene_names);         
        t2_expr_idx(~t2_expr_mask) = [];
        t2_targets(~t2_expr_mask) = [];
        
        t2_expr = Dream_expr(t2_expr_idx, t1_celllines);
        t1_expr = Dream_expr(t1_expr_idx, t2_celllines);

        
%         [~, SL_matches1] = ismember(t1_targets, SL_nodes);
%         t1_expr(~SL_matches1, :) = [];
%         SL_matches1(~SL_matches1) = [];
%         if(isempty(SL_matches1))
%             continue;
%         end
%         
%         [~, SL_matches2] = ismember(t2_targets, SL_nodes);
%         t2_expr(~SL_matches2, :) = [];
%         SL_matches2(~SL_matches2) = [];
%         if(isempty(SL_matches2))
%             continue;
%         end    

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
                
            p2(idx2) = ranksum(t2_Effect(active), t2_Effect(~active), 'tail', 'left');
        end
        p2_sum = Edgington(p2);
        
        pval = min(p1_sum, p2_sum);
%         pval2 = Edgington(x);
%         y = [sum(x), numel(x), pval, pval2]
%         pval = Edgington([p1; p2]);
                
        Vert_Syn_Confirm_pval(i, 1) = pval;        
    end
%%
    Vert_Syn_threshold = 0;

    Combined_SL_Syn = -log10(Vert_Syn_Confirm_pval) .* double(Vert_syn > Vert_Syn_threshold);
    Combined_SL_Syn = 100*(Combined_SL_Syn - min(Combined_SL_Syn)) ./ (max(Combined_SL_Syn) - min(Combined_SL_Syn));
        
%% Expression of functional neighborhood of each drug
%     Normalized_expr = zscore(Dream_expr, 0, 2); % DOES NOT WORK!
    Normalized_expr = Dream_expr;
    Mean_target_expr = zeros(size(annotations.drugs, 1), size(annotations.cellLines, 1));
    for i = 1:size(annotations.drugs, 1)
        targets = interactome.vertex_names(Drug_target_neighbors{i});

        [expr_mask, expr_idx] = ismember(targets, Dream_gene_names);     
        expr_idx(~expr_mask) = [];
        for j = 1:size(annotations.cellLines, 1)
            Mean_target_expr(i, j) = mean(Normalized_expr(expr_idx, j));
        end
    end
    
    Mean_target_expr_MidRange = abs(Mean_target_expr - 5.5);
%     Mean_target_expr_MidRange = (max(Mean_target_expr_MidRange(:)) - Mean_target_expr_MidRange) ./ (max(Mean_target_expr_MidRange(:)) - min(Mean_target_expr_MidRange(:)));
%% Validate
    DE = cell2mat(Mono.EMax);    
    DE(isempty(DE) | isnan(DE) | isinf(DE)) = 0; 
    DE = normalize(DE, 'dim', 1); 
    
    clear stats
    X = Pair_synergy;
    X(isinf(X)) = 0;
    for i = 1:size(Pairs, 1)
        [~, jj, vv] = find(X(i, :));
        stats(i, 1:3) = [nnz(vv > 20), nnz(vv> 30), nnz(vv>40)];
    end
    XXX = [Vert_syn, Horiz_syn, Combined_SL_Syn, stats];
    [Xc, Xp] = corr(XXX, 'type', 'Pearson')
    XT = mat2dataset(XXX);

    U = Combined_SL_Syn(Syn_pair_id);% .* max(Mean_target_expr_MidRange(sub2ind(size(Mean_target_expr), Syn_drug1, Syn_CL_id)), Mean_target_expr_MidRange(sub2ind(size(Mean_target_expr), Syn_drug2, Syn_CL_id)));
    U(isnan(U)) = 0;
    [CC, pp] = corr([U .* max(Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug1, Syn_CL_id)), Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug2, Syn_CL_id))), Syn_vv, Syn_vv > 20, Syn_vv > 30, Syn_vv > 40], 'type', 'Pearson')    
    
    
    Results = [annotations.drugs.ChallengeName(Syn_drug1), annotations.drugs.Class(Syn_drug1), num2cell(Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug1, Syn_CL_id))), num2cell(Mean_target_expr_MidRange(sub2ind(size(Mean_target_expr), Syn_drug1, Syn_CL_id))), annotations.drugs.ChallengeName(Syn_drug2), annotations.drugs.Class(Syn_drug2), num2cell(Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug2, Syn_CL_id))), num2cell(Mean_target_expr_MidRange(sub2ind(size(Mean_target_expr), Syn_drug2, Syn_CL_id))), num2cell(Combined_SL_Syn(Syn_pair_id)), num2cell(Syn_vv)];
    
    Results(Syn_vv < 20, :) = [];
    dlmcell('SL_results.tsv', Results);

%% Extract features for DT construction    
% 
%     Features_data = [Syn_vv > 30, Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug1, Syn_CL_id)), Mean_target_expr(sub2ind(size(Mean_target_expr), Syn_drug1, Syn_CL_id)), Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug2, Syn_CL_id)), Mean_target_expr(sub2ind(size(Mean_target_expr), Syn_drug2, Syn_CL_id)), Combined_SL_Syn(Syn_pair_id)];    
%     Features_names = {'DrugEffect1', 'NExpr1', 'DrugEffect2', 'NExpr2', 'Vert_Syn'};
%     Features_types = repmat({'continuous'}, numel(Features_names), 1);
% 
%     factor = 1;
%     pos_labels = find(Syn_labels == 1);
%     neg_labels = find(Syn_labels == 0);
% 
%     for i = 1:20
%         fd_names = fopen(sprintf('Dream_Sample%d.names', i), 'w');
%         fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
%         for k = 1:numel(Features_types)
%             fprintf(fd_names, '%s:\t%s.\n', Features_names{k}, Features_types{k});
%         end    
%         fclose(fd_names);
%         
%         neg_samples = randsample(neg_labels, numel(pos_labels)*factor);
%         samples = [pos_labels; neg_samples];
%         Sampled_Features = Features_data(samples, :);    
%         dlmcell(sprintf('Dream_Sample%d.data', i), num2cell(Sampled_Features), 'delimiter', ',');    
%     end    

%%
    D2D = Construct_D2D(annotations, interactome);
    [C2C, NodeWeights] = Construct_C2C(annotations, interactome, 'expression_only', true);
    
       
%     X = cell2mat(fillinBlanks(Mono.EMax, C2C, D2D));    
    X = cell2mat(Mono.EMax);    
    X(isempty(X) | isnan(X) | isinf(X)) = 0; 
    X = normalize(X, 'dim', 1);    

    V = -log10(Vert_Syn_Confirm_pval);
%     V(V<50) = 0;
    T = Modified_zscore(Mean_target_expr')';


%     X(X < 1/size(X, 1)) = 0;
    
    Z = zeros(size(Pairs, 1), size(annotations.cellLines, 1));
    for i = 1:size(Pairs, 1)
        for j = 1:size(annotations.cellLines, 1)
            Z(i, j) = V(i) * max([X(Pairs(i, 1), j), X(Pairs(i, 2), j)], [], 2);
        end
    end    
%     Z = Modified_zscore(Z')';
%     Z(isnan(Z)) = 0;
    outPath = fullfile('output', 'predictions', 'leaderBoard');    
    exportResults( annotations, Pairs, Pair_names, Z, outPath, 'synergy_threshold', 0.025);
    
%     ZZ = Z;
%     ZZ(ZZ < 0.055) = 0;
    %% TODO
    % Drug Effect: Currently EMax, but can be AUC, H, or IC50
    % Expression threshod: currently 3.9, but can vary
    % Confidence threshold on SL: currently 0.6, but can vary between 0-1
    % Threshold on Vert_Syn
    