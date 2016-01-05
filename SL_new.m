    clear
    addpath(genpath('code'));
    warning('off','all');    

    annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
    annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
    annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);

    [~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
    annotations.cellLines = annotations.cellLines(CL_perm, :);

    if(~exist('datasets.mat', 'file'))
        %% Read expression
        [expr_table, cellLine_names, Dream_gene_names] = my_tblread('input/Dream/molecular/gex.csv', ',');
        cellLine_names = cellfun(@(x) x(2:end-1), cellLine_names, 'UniformOutput', false);
        Dream_gene_names = cellfun(@(x) x(2:end-1), Dream_gene_names, 'UniformOutput', false);

        Dream_expr = zeros(size(expr_table, 1), size(annotations.cellLines, 1));    
        [~, selected_cols] = ismember(cellLine_names, annotations.cellLines.Sanger_Name);

        Dream_expr(:, selected_cols) = expr_table; %zscore(expr_table, 0, 2);

        drug_target_expr_indices = arrayfun(@(x) find(ismember(Dream_gene_names, annotations.drugs.Target{x})), 1:size(annotations.drugs, 1), 'UniformOutput', false);

        %% Read Monotherapy data
    %     fname = 'input/Dream/synergy/ch2_leaderBoard_monoTherapy.csv';
        fname = 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv';
    %     [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' ); %TODO: How to optimally combine replicates?

        [ Mono ] = read_singleMonoTherapy( annotations, fname );
        [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );

    %%
        synergy_threshold = 30;
        X = Pair_synergy;
        X(isinf(X)) = 0;
        [Syn_pair_id, Syn_CL_id, Syn_vv] = find(X);
        Syn_drug1 = Pairs(Syn_pair_id, 1);
        Syn_drug2 = Pairs(Syn_pair_id, 2);
        ck = Syn_CL_id;
        Syn_labels = synergy_threshold <= Syn_vv;    

    %% Read physical interactome
        [ interactome ] = readNetwork();

        drug_target_net_indices = arrayfun(@(x) find(ismember(interactome.vertex_names, annotations.drugs.Target{x})), 1:size(annotations.drugs, 1), 'UniformOutput', false);

    %% Read GO annotations and map them to the interactome
        fd = fopen('input/GO_subset/GO_annotations.txt', 'r');
        C = textscan(fd, '%[^\t] GO:%d %s');
        fclose(fd);

        [vertex_mask, vertex_id] = ismember(C{1}, interactome.vertex_genes);
        vertex_id = vertex_id(vertex_mask);
        GO_id = C{2}(vertex_mask);
        GO_annotations = cell(numel(interactome.vertex_genes), 1);
        for i = 1:numel(GO_id)
            GO_annotations{vertex_id(i)} = [GO_annotations{vertex_id(i)}; GO_id(i)];
        end
        total_GO = unique(GO_id);
        interactome.vertex_annotations = GO_annotations;

        pop_size = numel(total_GO);
        SemSim = zeros(numel(interactome.vertex_genes));
        for v = 1:numel(interactome.vertex_genes)
            fprintf('%d\n', v);
            for u = v+1:numel(interactome.vertex_genes)
                success = numel(intersect(interactome.vertex_annotations{u}, interactome.vertex_annotations{v}));
                SemSim(u, v) = -log10(hygecdf(success, pop_size, numel(interactome.vertex_annotations{u}), numel(interactome.vertex_annotations{v}), 'upper'));
            end        
        end
        SemSim = max(SemSim, SemSim');    
        SemSim(isinf(SemSim)) = 0;

    %% Expand neighborhood of interactome vertices
        interactome.Neighbors = arrayfun(@(v) find(interactome.A(v, :)), 1:numel(interactome.vertex_genes), 'UniformOutput', false)';    

        interactome.Pruned_neighbors = arrayfun(@(v) union(v, interactome.Neighbors{v}(arrayfun(@(u) SemSim(u,v) > 10 , interactome.Neighbors{v}))), 1:numel(interactome.vertex_genes), 'UniformOutput', false)';




        Drug_target_neighbors = cellfun(@(d) unique([interactome.Pruned_neighbors{d}]), drug_target_net_indices, 'UniformOutput', false)';
    
    else
        load datasets
    end
    
%%
    Prctl = cell2mat(arrayfun(@(c) prctile(Dream_expr(:, c), 1:1:100), 1:size(Dream_expr, 2), 'UniformOutput', false)');
    expr_thresholds = median(Prctl);
    Dream_expr_digitized = zeros(size(Dream_expr));
    Dream_expr_digitized(Dream_expr <= expr_thresholds(1)) = 1;
    for i = 1:numel(expr_thresholds)-1
        Dream_expr_digitized(expr_thresholds(i) < Dream_expr & Dream_expr <= expr_thresholds(i+1)) = i+1;        
    end

%%    
    Drug_Effect = cell2mat(Mono.EMax);
    Drug_Effect(isnan(Drug_Effect)) = 100;
    Drug_Effect = 100 - Drug_Effect; % effectiveness of drug in cell-lines

%     Drug_Effect = cell2mat(Mono.Drug_sensitivity);
%     Drug_Effect(isnan(Drug_Effect)) = 0;

    expression_threshold = 5;
    
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
    [CC, pp] = corr([U>20, Syn_vv > 30], 'type', 'Pearson')    

    
    
    Results = [annotations.drugs.ChallengeName(Syn_drug1), annotations.drugs.MoA(Syn_drug1), num2cell(Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug1, Syn_CL_id))), num2cell(Mean_target_expr_MidRange(sub2ind(size(Mean_target_expr), Syn_drug1, Syn_CL_id))), annotations.drugs.ChallengeName(Syn_drug2), annotations.drugs.MoA(Syn_drug2), num2cell(Drug_Effect(sub2ind(size(Drug_Effect), Syn_drug2, Syn_CL_id))), num2cell(Mean_target_expr_MidRange(sub2ind(size(Mean_target_expr), Syn_drug2, Syn_CL_id))), num2cell(Combined_SL_Syn(Syn_pair_id)), num2cell(Syn_vv)];
    
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
    V = Combined_SL_Syn;
    V(V<30) = 0;
    
    Z = zeros(size(Pairs, 1), size(annotations.cellLines, 1));
    for i = 1:size(Pairs, 1)
        for j = 1:size(annotations.cellLines, 1)
            Z(i, j) = V(i);% * (max(Mean_target_expr_MidRange(Pairs(i, 1), j), Mean_target_expr_MidRange(Pairs(i, 2), j)) < 0.75);
        end
    end
%     Z = Modified_zscore(Z')';
%     Z(isnan(Z)) = 0;
    outPath = fullfile('output', 'predictions', 'leaderBoard');    
    exportResults( annotations, Pairs, Pair_names, Z, outPath, 'synergy_threshold', 0);
    %% TODO
    % Drug Effect: Currently EMax, but can be AUC, H, or IC50
    % Expression threshod: currently 3.9, but can vary
    % Confidence threshold on SL: currently 0.6, but can vary between 0-1
    % Threshold on Vert_Syn
    