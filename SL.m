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
    [Syn_pair_id, Syn_CL_id, Syn_vv] = find(X);
    Syn_drug1 = Pairs(Syn_pair_id, 1);
    Syn_drug2 = Pairs(Syn_pair_id, 2);
    ck = Syn_CL_id;
    Syn_labels = synergy_threshold <= Syn_vv;

%%
    Syn_vals = sort(Syn_vv(Syn_vv < prctile(Syn_vv, 95) & Syn_vv > prctile(Syn_vv, 5)));
    fitdist(Syn_vals,'tlocationscale')
    
    prctile(Syn_vv, 5:5:95)
%%
    [ interactome ] = readNetwork();

%% Drug-Drug synthetic lethal feature (parallel pathways)
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

        [~, direct_SL_matches1{i}] = ismember(t1, SL_nodes);
        direct_SL_matches1{i}(~direct_SL_matches1{i}) = [];
        
        N1 = [];
        t1_subnet = interactome.A(t1_idx, :);
        for j = 1:numel(t1_idx)
            N1 = union(N1, find(t1_subnet(j, :)));
        end
        N1 = union(t1_idx, N1);
        
        if(numel(N1) >= 2)            
            N1_density(i) = 100*nnz(interactome.A(N1, N1))/2 / nchoosek(numel(N1), 2);
        else
            N1_density(i) = 0;
        end
        N1_names = interactome.vertex_names(N1);
        [~, SL_matches1] = ismember(N1_names, SL_nodes);
        SL_matches1(~SL_matches1) = [];


        t2 = annotations.drugs.Target{Pairs(i, 2)};
        [~, t2_idx] = ismember(t2, interactome.vertex_names); 
        t2_idx(~t2_idx) = [];
        t2(~t2_idx) = [];
        
        [~, direct_SL_matches2{i}] = ismember(t2, SL_nodes);
        direct_SL_matches2{i}(~direct_SL_matches2{i}) = [];
        

        N2 = [];
        t2_subnet = interactome.A(t2_idx, :);        
        for j = 2:numel(t2_idx)
            N2 = union(N2, find(t2_subnet(j, :)));
        end
        N2 = union(t2_idx, N2);
        if(numel(N2) >=2)
            N2_density(i) = 100*nnz(interactome.A(N2, N2))/2 / nchoosek(numel(N2), 2);
        else
            N2_density(i) = 0;
        end
        
        N2_names = interactome.vertex_names(N2);
        [~, SL_matches2] = ismember(N2_names, SL_nodes);
        SL_matches2(~SL_matches2) = [];


        
        success_size(i) = nnz(SL_net(SL_matches1, SL_matches2));
        sample_size = numel(SL_matches1)*numel(SL_matches2);

        if(~isempty(direct_SL_matches1{i}) && ~isempty(direct_SL_matches1{i}))
            direct_Syn(i) = nnz(SL_net(direct_SL_matches1{i}, direct_SL_matches2{i}));
        else
            direct_Syn(i) = 0;
        end
        Syn(i) = -log10(min(1, hygecdf(success_size(i), population_size, total_success, sample_size, 'upper') * size(Pairs, 1)));
    end
    Syn_norm = Syn;    
    Syn_norm(isinf(Syn_norm)) = 1.5*max(nonzeros(Syn_norm(~isinf(Syn_norm))));
    Syn_norm = Syn_norm ./ nanstd(Syn_norm);

    XXX = [Syn_vv, direct_Syn(Syn_pair_id)', N1_density(Syn_drug1)', N2_density(Syn_drug2)', Syn_norm(Syn_pair_id), [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)]', Syn_drug1, Syn_CL_id)}]' , [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], Syn_drug2, Syn_CL_id)}]']
%     XXX = zscore(XXX)
    XXX_table = mat2dataset(XXX);
    
    Syn_norm = Syn;    
    Syn_norm(isinf(Syn_norm)) = 1.5*max(nonzeros(Syn_norm(~isinf(Syn_norm))));
    Syn_norm = Syn_norm ./ nanstd(Syn_norm);
    features.Drug.data = [features.Drug.data, num2cell(Syn_norm(Syn_pair_id))];
    features.Drug.name = [features.Drug.name; 'Syn'];
    features.Drug.type = [features.Drug.type; 'continuous'];

