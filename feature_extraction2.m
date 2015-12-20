    clear
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
    [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' ); %TODO: How to optimally combine replicates?

%     [ Mono ] = read_singleMonoTherapy( annotations, fname );
    [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );

%%
    synergy_threshold = 30;
    X = Pair_synergy;
    X(isinf(X)) = 0;
    [Syn_ii, Syn_jj, Syn_vv] = find(X);
    di = Pairs(Syn_ii, 1);
    dj = Pairs(Syn_ii, 2);
    ck = Syn_jj;
    Syn_labels = synergy_threshold <= Syn_vv;

%% Drug-Drug primary features
    [ interactome ] = readNetwork();

    fprintf('Computing Synergy for drug pairs based on SL interactions\n');
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
        [~, t1_idx] = ismember(t1, interactome.vertex_names); 
        t1_idx(~t1_idx) = [];
        t1(~t1_idx) = [];
        t1_subnet = interactome.A(t1_idx, :);
        N1 = [];
        for j = 1:numel(t1_idx)
            N1 = union(N1, find(t1_subnet(j, :)));
        end
        N1 = union(t1_idx, N1);
        N1_names = interactome.vertex_names(N1);
        [~, SL_matches1] = ismember(N1_names, SL_nodes);
        SL_matches1(~SL_matches1) = [];


        t2 = annotations.drugs.Target{Pairs(i, 2)};
        [~, t2_idx] = ismember(t2, interactome.vertex_names); 
        t2_idx(~t2_idx) = [];
        t2(~t2_idx) = [];
        t2_subnet = interactome.A(t2_idx, :);
        N2 = [];
        for j = 2:numel(t2_idx)
            N2 = union(N2, find(t2_subnet(j, :)));
        end
        N2 = union(t2_idx, N2);
        N2_names = interactome.vertex_names(N2);
        [~, SL_matches2] = ismember(N2_names, SL_nodes);
        SL_matches2(~SL_matches2) = [];

        success_size = nnz(SL_net(SL_matches1, SL_matches2));
        sample_size = numel(SL_matches1)*numel(SL_matches2);


        Syn(i) = -log10(min(1, hygecdf(success_size, population_size, total_success, sample_size, 'upper') * size(Pairs, 1)));
    end

    Syn_norm = Syn < 150;    
%     Syn_norm(isinf(Syn_norm)) = 10*max(nonzeros(Syn_norm(~isinf(Syn_norm))));
%     Syn_norm = Syn_norm ./ nanstd(Syn_norm);

    Features_data_mat(:, 1) = Syn_norm(Syn_ii);
	Features_data_mat(:, 2) = max([Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Features_data_mat(:, 3) = min([Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Features_data_mat(:, 4) = max([Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Features_data_mat(:, 5) = min([Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
%     Features_data_mat = zscore(Features_data_mat);
    Features_data = num2cell(Features_data_mat);
    
    Features_name = {'SL', 'MaxDrug_sensitivity', 'MinDrug_sensitivity', 'MaxMaxC', 'MinMaxC'};
    Features_type = repmat({'continuous'}, 5, 1);

    %%
    factor = 1;
    pos_labels = find(Syn_labels == 1);
    neg_labels = find(Syn_labels == 0);

    for i = 1:20
        fd_names = fopen(sprintf('Dream_Sample%d.names', i), 'w');
        fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
        for k = 1:numel(Features_name)
            fprintf(fd_names, '%s:\t%s.\n', Features_name{k}, Features_type{k});
        end    
        fclose(fd_names);
        
        neg_samples = randsample(neg_labels, numel(pos_labels)*factor);
        samples = [pos_labels; neg_samples];
        Sampled_Features = [num2cell(Syn_labels(samples, 1)), Features_data(samples, :)];    
        dlmcell(sprintf('Dream_Sample%d.data', i), Sampled_Features, 'delimiter', ',');    
    end
