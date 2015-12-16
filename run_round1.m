    clear
    addpath(genpath('code'));
    warning('off','all');    

    annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
    annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
    annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);

    [~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
    annotations.cellLines = annotations.cellLines(CL_perm, :);
    
%% Load signaling interactome from ACSN together with 55 functional classes relevant to cancer
    % TODO: Edge directions
    [ ACSN ] = import_ACSN();


%% Read Monotherapy data and impute missing values based on Dual Layer method
%     fname = 'input/Dream/synergy/ch2_leaderBoard_monoTherapy.csv';
    fname = 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv';
    [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' );

%     [ Mono ] = read_singleMonoTherapy( annotations, fname );
    [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );


%% Drug-Drug Similarity Network
    D2D = Construct_D2D(annotations, ACSN);

    
%% Cellline-Celline Similarity Network
    % TODO: Which networks to use? Should we also use co-methyl and co-mut? How
    % to optimally combine using Mashup, GeneMANIA, or SNF?
    [C2C, NodeWeights] = Construct_C2C(annotations, ACSN, 'expression_only', true);
    
%% Fill in the gaps
%     fprintf('Imputing missing elements in Monotherapy data\n');
%     Mono.IC50 = fillinBlanks(Mono.IC50, C2C, D2D);
%     Mono.H = fillinBlanks(Mono.H, C2C, D2D);
%     Mono.EMax = fillinBlanks(Mono.EMax, C2C, D2D);
%     Mono.Max_C = fillinBlanks(Mono.Max_C, C2C, D2D);
%     Mono.Drug_sensitivity = fillinBlanks(Mono.Drug_sensitivity, C2C, D2D);
    
%%
    fprintf('Computing Synergy based on SL interactions\n');
    fd = fopen('/home/shahin/Dropbox/Dream/experiment/refs/GeneticInteractions/SynLethDB_human.txt', 'r');
    C = textscan(fd, '%s %s');
    fclose(fd);
    SL_nodes = union(C{1}, C{2});
    [~, ii] = ismember(C{1}, SL_nodes);
    [~, jj] = ismember(C{2}, SL_nodes);
    SL_net = sparse(ii, jj, 1, numel(SL_nodes), numel(SL_nodes));
    SL_net = max(SL_net, SL_net');

    population_size = nchoosek(numel(SL_nodes), 2);
    total_success = nnz(SL_net)/2;

    Syn = zeros(size(Pairs, 1), 1);
    for i = 1:size(Pairs, 1)
        t1 = annotations.drugs.Target{Pairs(i, 1)};
        [~, t1_idx] = ismember(t1, ACSN.vertex_names); 
        t1_idx(~t1_idx) = [];
        t1(~t1_idx) = [];
        t1_subnet = ACSN.A(t1_idx, :);
        N1 = [];
        for j = 1:numel(t1_idx)
            N1 = union(N1, find(t1_subnet(j, :)));
        end
        N1 = t1_idx; %union(t1_idx, N1);
        N1_names = ACSN.vertex_names(N1);
        [~, SL_matches1] = ismember(N1_names, SL_nodes);
        SL_matches1(~SL_matches1) = [];


        t2 = annotations.drugs.Target{Pairs(i, 2)};
        [~, t2_idx] = ismember(t2, ACSN.vertex_names); 
        t2_idx(~t2_idx) = [];
        t2(~t2_idx) = [];
        t2_subnet = ACSN.A(t2_idx, :);
        N2 = [];
        for j = 2:numel(t2_idx)
            N2 = union(N2, find(t2_subnet(j, :)));
        end
        N2 = t2_idx; %union(t2_idx, N2);
        N2_names = ACSN.vertex_names(N2);
        [~, SL_matches2] = ismember(N2_names, SL_nodes);
        SL_matches2(~SL_matches2) = [];

        success_size = nnz(SL_net(SL_matches1, SL_matches2));
        sample_size = numel(SL_matches1)*numel(SL_matches2);


        Syn(i) = success_size; %-log10(min(1, hygecdf(success_size, population_size, total_success, sample_size, 'upper') * size(Pairs, 1)));
    end


%%
    Syn_norm = Syn;    
    Syn_norm(isinf(Syn_norm)) = 1.5*max(nonzeros(Syn_norm(~isinf(Syn_norm))));
    Syn_norm = Syn_norm ./ nanstd(Syn_norm);

    Syn_w = ones(size(Syn));
    Syn_w(Syn > 10) = 16;
    
    Z = Pair_synergy;
    Z(Z < 0) = 0;
    [ii, jj, vv] = find(Z);
    Z_ind = sub2ind(size(Z), ii, jj); % for validating Y

    clear ZZ
    for i = 1:size(Z, 1)
        row = nonzeros(Z(i, :));
        val = min(row);
        if(~isempty(val))
            ZZ(i, 1) = val;    
        else
            ZZ(i, 1) = 0;
        end

        val = median(row);
        if(~isempty(val))
            ZZ(i, 2) = val;    
        else
            ZZ(i, 2) = 0;
        end
        
        val = max(row);
        if(~isempty(val))
            ZZ(i, 3) = val;    
        else
            ZZ(i, 3) = 0;
        end
        
        val = nnz(row > 30);
        if(~isempty(val))
            ZZ(i, 4) = val;    
        else
            ZZ(i, 4) = 0;
        end        
        
        ZZ(i, 5) = nanmedian(nonzeros(cell2mat(Mono.IC50(Pairs(i, 1), :))));    
        ZZ(i, 6) = nanmedian(nonzeros(cell2mat(Mono.IC50(Pairs(i, 2), :))));
        ZZ(i, 7) = Syn(i);
        
    %     ZZ(i) = max(Z(i, :));    
    end    
    ZZ(:, 5) = - zscore(ZZ(:, 5));
    ZZ(:, 6) = - zscore(ZZ(:, 6));
    
    ZZ(isnan(ZZ) | isinf(ZZ)) = 0;
    

    [c, p] = corr(Syn_norm, ZZ)

    %%
    [ transcriptional_gene_signature, transc_class_gene_idx ] = computeTransSig( annotations, ACSN );
    
    for i = 1:size(Pairs, 1)
        D1 = Pairs(i, 1);
        D2 = Pairs(i, 2);
        if( isempty(transcriptional_gene_signature{D1, 2}) || isempty(transcriptional_gene_signature{D2, 2}) )
            T(i, 1) = 1;
        else
            T(i, 1)= 1 + corr(transcriptional_gene_signature{D1, 2}, transcriptional_gene_signature{D2, 2});
        end
    end
    
%%

    X = cell2mat(Mono.Drug_sensitivity);    
    X(isempty(X) | isnan(X) | isinf(X)) = 0;   
    X = normalize(X, 'dim', 1);
    
    clear Y;
    for i = 1:size(Pairs, 1)
        for j = 1:size(annotations.cellLines, 1)
            Y(i, j) = (X(Pairs(i, 1), j) * X(Pairs(i, 2), j));
%             Y(i, j) = max([X(Pairs(i, 1), j), X(Pairs(i, 2), j)], [], 2);
        end
    end
    
%     Y = bsxfun(@times, Y, T);
    Y = Y ./ nanstd(nonzeros(Y));    
    Y = bsxfun(@times, Y, Syn_w);
%     Y = Y ./ max(Y(:));
    
%     [ii, jj, vv] = find(Z);
%     ind = sub2ind(size(Z), ii, jj);
%     corr(vv, Y(ind))
%     [ii, jj, vv] = find(Z);
%     vv_rank = tiedrank(vv);
%     Z_rank = full(sparse(ii, jj, vv_rank, size(Z, 1), size(Z, 2)));
%     Y = reshape(tiedrank(Y(:)), size(Y, 1), size(Y, 2));

    
    [cc, pval] = corr(Z(Z_ind), Y(Z_ind), 'type', 'Pearson')

    outPath = fullfile('output', 'predictions', 'leaderBoard');    
    exportResults( annotations, Pairs, Pair_names, Y, outPath, 'synergy_threshold', prctile(nonzeros(Y), 80));


