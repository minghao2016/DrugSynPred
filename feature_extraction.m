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
    
%% Celline features
    all_tissues = unique(annotations.cellLines.Tissue__General_);    
    [~, CL_tissue_labels] = ismember(annotations.cellLines.Disease_Area, all_tissues);
%     features.CL(1:numel(ck), 1) = Disease_Area(ck); % tissue of origin
    features.CL.data(1:numel(ck), 1) = annotations.cellLines.Disease_Area(ck);
    features.CL.name = {'Tissue of origin'};
    features.CL.type = {strjoin(unique(annotations.cellLines.Disease_Area)', ', ')};
    
%% Drug features   
    MW_data(:, 1) = abs(annotations.drugs.MW(di) - annotations.drugs.MW(dj));
    MW_data(:, 2) = max([annotations.drugs.MW(di) , annotations.drugs.MW(dj)], [], 2);
    MW_data(:, 3) = min([annotations.drugs.MW(di) , annotations.drugs.MW(dj)], [], 2);
    MW_data(:, 4) = annotations.drugs.MW(di) + annotations.drugs.MW(dj);
    MW_data(:, 5) = annotations.drugs.MW(di) .* annotations.drugs.MW(dj);
    MW_name = {'DeltaMW', 'MaxMW', 'MinMW', 'SumMW', 'ProdMW'};
    
    HBA_data(:, 1) = abs(annotations.drugs.HBA(di) - annotations.drugs.HBA(dj));
    HBA_data(:, 2) = max([annotations.drugs.HBA(di) , annotations.drugs.HBA(dj)], [], 2);
    HBA_data(:, 3) = min([annotations.drugs.HBA(di) , annotations.drugs.HBA(dj)], [], 2);
    HBA_data(:, 4) = annotations.drugs.HBA(di) + annotations.drugs.HBA(dj);
    HBA_data(:, 5) = annotations.drugs.HBA(di) .* annotations.drugs.HBA(dj);
    HBA_name = {'DeltaHBA', 'MaxHBA', 'MinHBA', 'SumHBA', 'ProdHBA'};

    HBD_data(:, 1) = abs(annotations.drugs.HBD(di) - annotations.drugs.HBD(dj));
    HBD_data(:, 2) = max([annotations.drugs.HBD(di) , annotations.drugs.HBD(dj)], [], 2);
    HBD_data(:, 3) = min([annotations.drugs.HBD(di) , annotations.drugs.HBD(dj)], [], 2);
    HBD_data(:, 4) = annotations.drugs.HBD(di) + annotations.drugs.HBD(dj);
    HBD_data(:, 5) = annotations.drugs.HBD(di) .* annotations.drugs.HBD(dj);
    HBD_name = {'DeltaHBD', 'MaxHBD', 'MinHBD', 'SumHBD', 'ProdHBD'};
    
    cLogP_data(:, 1) = abs(annotations.drugs.cLogP(di) - annotations.drugs.cLogP(dj));
    cLogP_data(:, 2) = max([annotations.drugs.cLogP(di) , annotations.drugs.cLogP(dj)], [], 2);
    cLogP_data(:, 3) = min([annotations.drugs.cLogP(di) , annotations.drugs.cLogP(dj)], [], 2);
    cLogP_data(:, 4) = annotations.drugs.cLogP(di) + annotations.drugs.cLogP(dj);
    cLogP_data(:, 5) = annotations.drugs.cLogP(di) .* annotations.drugs.cLogP(dj);
    
    Lipinski_data(:, 1) = abs(annotations.drugs.Lipinski(di) - annotations.drugs.Lipinski(dj));
    Lipinski_data(:, 2) = max([annotations.drugs.Lipinski(di) , annotations.drugs.Lipinski(dj)], [], 2);
    Lipinski_data(:, 3) = min([annotations.drugs.Lipinski(di) , annotations.drugs.Lipinski(dj)], [], 2);
    Lipinski_data(:, 4) = annotations.drugs.Lipinski(di) + annotations.drugs.Lipinski(dj);
    Lipinski_data(:, 5) = annotations.drugs.Lipinski(di) .* annotations.drugs.Lipinski(dj);
	Lipinski_name = {'DeltaLipinski', 'MaxLipinski', 'MinLipinski', 'SumLipinski', 'ProdLipinski'};
	cLogP_name = {'DeltacLogP', 'MaxcLogP', 'MincLogP', 'SumcLogP', 'ProdcLogP'};
    
    features.Drug.data = num2cell([MW_data, HBA_data, HBD_data, cLogP_data, Lipinski_data]);
    features.Drug.name = [MW_name, HBA_name, HBD_name, cLogP_name, Lipinski_name]';
    features.Drug.type = repmat({'continuous'}, 25, 1);
    
%%  Drug-CellLine features

	IC50_data(:, 1) = abs([Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' - [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]');
	IC50_data(:, 2) = max([Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	IC50_data(:, 3) = min([Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	IC50_data(:, 4) = [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' + [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	IC50_data(:, 5) = [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' .* [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	IC50_name = {'DeltaIC50', 'MaxIC50', 'MinIC50', 'SumIC50', 'ProdIC50'};

	H_data(:, 1) = abs([Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' - [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]');
	H_data(:, 2) = max([Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	H_data(:, 3) = min([Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	H_data(:, 4) = [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' + [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	H_data(:, 5) = [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' .* [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	H_name = {'DeltaH', 'MaxH', 'MinH', 'SumH', 'ProdH'};

	EMax_data(:, 1) = abs([Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' - [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]');
	EMax_data(:, 2) = max([Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	EMax_data(:, 3) = min([Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	EMax_data(:, 4) = [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' + [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	EMax_data(:, 5) = [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' .* [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	EMax_name = {'DeltaEMax', 'MaxEMax', 'MinEMax', 'SumEMax', 'ProdEMax'};

	Max_C_data(:, 1) = abs([Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' - [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]');
	Max_C_data(:, 2) = max([Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Max_C_data(:, 3) = min([Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Max_C_data(:, 4) = [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' + [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	Max_C_data(:, 5) = [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' .* [Mono.Max_C{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	Max_C_name = {'DeltaMax_C', 'MaxMax_C', 'MinMax_C', 'SumMax_C', 'ProdMax_C'};

	Drug_sensitivity_data(:, 1) = abs([Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' - [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]');
	Drug_sensitivity_data(:, 2) = max([Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Drug_sensitivity_data(:, 3) = min([Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}] , [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}])';
	Drug_sensitivity_data(:, 4) = [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' + [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	Drug_sensitivity_data(:, 5) = [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]' .* [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
	Drug_sensitivity_name = {'DeltaDrug_sensitivity', 'MaxDrug_sensitivity', 'MinDrug_sensitivity', 'SumDrug_sensitivity', 'ProdDrug_sensitivity'};

    features.DrugCL.data = num2cell([IC50_data, H_data, EMax_data, Max_C_data, Drug_sensitivity_data]);
    features.DrugCL.name = [IC50_name, H_name, EMax_name, Max_C_name, Drug_sensitivity_name]';
    features.DrugCL.type = repmat({'continuous'}, 25, 1);
%     
% %% Drug-Drug primary features
    [ interactome ] = readNetwork();
% 
%     [ D2D, coreD2D, targetD2D ] = Construct_D2D(annotations, interactome); % 1) How to optimally combine core and target D2D?, 2) Make sure Stitch core network is correct, 3) Target extraction using Stich, WinDTome, ...
%     features.Drug1Drug2(1:numel(di), 1) = full(coreD2D(sub2ind(size(coreD2D), di, dj)));
%     features.Drug1Drug2(1:numel(di), 2) = full(targetD2D(sub2ind(size(targetD2D), di, dj)));
% 
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
        t1_subnet = interactome.A(t1_idx, :);
        N1 = [];
        for j = 1:numel(t1_idx)
            N1 = union(N1, find(t1_subnet(j, :)));
        end
        N1 = t1_idx; %union(t1_idx, N1);
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
        N2 = t2_idx; %union(t2_idx, N2);
        N2_names = interactome.vertex_names(N2);
        [~, SL_matches2] = ismember(N2_names, SL_nodes);
        SL_matches2(~SL_matches2) = [];

        success_size = nnz(SL_net(SL_matches1, SL_matches2));
        sample_size = numel(SL_matches1)*numel(SL_matches2);


        Syn(i) = -log10(min(1, hygecdf(success_size, population_size, total_success, sample_size, 'upper') * size(Pairs, 1)));
    end

    Syn_norm = Syn;    
    Syn_norm(isinf(Syn_norm)) = 1.5*max(nonzeros(Syn_norm(~isinf(Syn_norm))));
    Syn_norm = Syn_norm ./ nanstd(Syn_norm);
    features.Drug.data = [features.Drug.data, num2cell(Syn_norm(Syn_ii))];
    features.Drug.name = [features.Drug.name; 'Syn'];
    features.Drug.type = [features.Drug.type; 'continuous'];
%     
%     features.Drug1Drug2(1:numel(di), 3) = Syn_norm(Syn_ii);
%     
% %% Drug-Drug-CL
%     % a) SL_sim only using expressed genes
%     
%     % b) Similarity of gene expression: 
%     % 1) Read LINC dataset
%     % 2) Compute Sim/Dissimilarity (different measures)
%     % 3) Project onto Census genes or ACSN functional classes before sim/dissim compitations
     
%% Train a predictor to combine features using Random Forest
%     Feature_data = [features.CL.data, features.Drug.data, features.DrugCL.data];
%     Feature_name = [features.CL.name; features.Drug.name; features.DrugCL.name];    
%     Feature_type = [features.CL.type; features.Drug.type; features.DrugCL.type];    
    Feature_data = [features.DrugCL.data];
    Feature_name = [features.DrugCL.name];    
    Feature_type = [features.DrugCL.type];    

    
%     fd_names = fopen(sprintf('Dream.names', i), 'w');
%     fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
%     for k = 1:numel(Feature_name)
%         fprintf(fd_names, '%s:\t%s.\n', Feature_name{k}, Feature_type{k});
%     end    
%     fclose(fd_names);
% 
%     dlmcell('Dream.data', [num2cell(Syn_labels), Feature_data], 'delimiter', ',');    
%     

    
    
    factor = 2;
    pos_labels = find(Syn_labels == 1);
    neg_labels = find(Syn_labels == 0);

    for i = 1:20
        fd_names = fopen(sprintf('Dream_Sample%d.names', i), 'w');
        fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
        for k = 1:numel(Feature_name)
            fprintf(fd_names, '%s:\t%s.\n', Feature_name{k}, Feature_type{k});
        end    
        fclose(fd_names);
        
        neg_samples = randsample(neg_labels, numel(pos_labels)*factor);
        samples = [pos_labels; neg_samples];
        Sampled_Features = [num2cell(Syn_labels(samples, 1)), Feature_data(samples, :)];    
        dlmcell(sprintf('Dream_Sample%d.data', i), Sampled_Features, 'delimiter', ',');    
    end
    
%%    
    Drug_Range = 1 + 2:11;
    DefaultTree = fitctree(Features(:, Drug_Range), Syn_labels, 'PredictorNames', Feature_labels(Drug_Range), 'MinLeaf', 50, 'PruneCriterion', 'error');
    view(DefaultTree,'Mode','Graph')
    
    DDC_Range = [24];
    DefaultTree = fitctree(Features(:, DDC_Range), Syn_labels, 'PredictorNames', Feature_labels(DDC_Range));
    view(DefaultTree,'Mode','Graph')
    
%%
    RF_ensemble = TreeBagger(500, Features, Syn_labels,'Method','Classification','OOBVarImp','On');
%     RF_ensemble.OOBPermutedVarDeltaError  
    
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    h = bar(RF_ensemble.OOBPermutedVarDeltaError);

    set(gca, 'XLim', [0, numel(RF_ensemble.OOBPermutedVarDeltaError) + 1]);
    set(gca,'FontSize', 16, 'FontWeight','bold'); 
    xticklabel_rotate(1:numel(Feature_labels), 45, Feature_labels, 'FontSize', 16, 'FontWeight','bold', 'Interpreter', 'none');        
    set(h(1),'FaceColor', [0.5 0.5 0.5]);
    print(fig, 'output/plots/feature_importance.eps', '-depsc2','-r300');
    
    
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    oobErrorBaggedEnsemble = oobError(RF_ensemble);
    plot(oobErrorBaggedEnsemble, 'LineWidth', 2)
    set(gca,'FontSize', 14, 'FontWeight','bold'); 
    xlabel( 'Number of grown trees','FontSize', 16, 'FontWeight','bold' );
    ylabel('Out-of-bag classification error','FontSize', 16, 'FontWeight','bold' );   
    
    print(fig, 'output/plots/out-of-bag_error.eps', '-depsc2','-r300');
    
%% Construct optimal Decision Tree and choose the right depth
    leafs = logspace(1,2,10); % Generate minimum leaf occupancies for classification trees from 10 to 100, spaced exponentially apart

    rng('default')
    N = numel(leafs);
    err = zeros(N,1);
    for n=1:N
        t = fitctree(Features, Syn_labels, 'CrossVal', 'On', 'MinLeaf', leafs(n));
        err(n) = kfoldLoss(t);
    end
    plot(leafs,err);
    xlabel('Min Leaf Size');
    ylabel('cross-validated error');


%%
DefaultTree = fitctree(Features, Syn_labels, 'PredictorNames', Feature_labels, 'MinLeaf', 40);
view(DefaultTree,'Mode','Graph')

OptimalTree = fitctree(Features, Syn_labels, 'minleaf', 20, 'PredictorNames', Feature_labels);
view(OptimalTree,'mode','graph')

resubOpt = resubLoss(OptimalTree);
lossOpt = kfoldLoss(crossval(OptimalTree));
resubDefault = resubLoss(DefaultTree);
lossDefault = kfoldLoss(crossval(DefaultTree));
resubOpt,resubDefault,lossOpt,lossDefault

%%
    OptimalTree = fitctree(Features, Syn_labels, 'MinLeaf', 30, 'PredictorNames', Feature_labels);
    view(OptimalTree, 'mode','graph')
