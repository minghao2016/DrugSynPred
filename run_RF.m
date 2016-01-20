    clear
    addpath(genpath('code'));
    warning('off','all');    

    annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
    annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
    annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);

    [~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
    annotations.cellLines = annotations.cellLines(CL_perm, :);

        %% Read Monotherapy data
%         fname = 'input/Dream/synergy/ch2_leaderBoard_monoTherapy.csv';
        fname = 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv';
        [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' ); %TODO: How to optimally combine replicates?

%         [ Mono ] = read_singleMonoTherapy( annotations, fname );
        [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );
%%
%     [ interactome ] = readNetwork();
    load('datasets', 'interactome');
    D2D = D2D_combined(annotations, interactome); % Construct_D2D(annotations, interactome);
%     [C2C, NodeWeights] =  Construct_C2C(annotations, interactome, 'expression_only', true);    
    C2C =  C2C_expr( annotations );
    Mono.EMax = fillinBlanks(Mono.EMax, C2C, D2D);

     

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

%% Drug-Drug primary features
    % Vertical (topology)
    [ VertSyn_Topo ] = computeVertSyn_topo( annotations, interactome, Pairs );
    features.DrugDrug.data = num2cell(VertSyn_Topo(Syn_ii));
    features.DrugDrug.name = {'VertSyn_Topo'};
    features.DrugDrug.type = {'continuous'};
    
    % Vertical (expression)
    [ VertSyn_Expr, VertSyn_Expr_Rest ] = computeVertSyn_expr( annotations, Pairs );
    features.DrugDrug.data = [features.DrugDrug.data, num2cell(VertSyn_Expr(Syn_ii))];
    features.DrugDrug.name = [features.DrugDrug.name; 'VertSyn_Expr'];
    features.DrugDrug.type = [features.DrugDrug.type; 'continuous'];
    
    % Horizental (SL)
    Drug_Effect = cell2mat(Mono.EMax);
    Drug_Effect(isnan(Drug_Effect)) = 100;
    Drug_Effect = 100 - Drug_Effect; % effectiveness of drug in cell-lines
    
    [ Horiz_compSL_Syn ] = computeHorSyn_DrugEffect( annotations, Pairs, Drug_Effect );  
    features.DrugDrug.data = [features.DrugDrug.data, num2cell(Horiz_compSL_Syn(Syn_ii))];
    features.DrugDrug.name = [features.DrugDrug.name; 'HorizSyn_Effect'];
    features.DrugDrug.type = [features.DrugDrug.type; 'continuous'];    
    
    %

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
    Feature_data = [features.DrugCL.data, features.DrugDrug.data];
    Feature_name = [features.DrugCL.name; features.DrugDrug.name];    
    Feature_type = [features.DrugCL.type; features.DrugDrug.type];  

    
%     Feature_data = [features.CL.data, features.Drug.data, features.DrugCL.data, features.DrugDrug.data];
%     Feature_name = [features.CL.name; features.Drug.name; features.DrugCL.name; features.DrugDrug.name];    
%     Feature_type = [features.CL.type; features.Drug.type; features.DrugCL.type; features.DrugDrug.type];    
%     Feature_data = [features.DrugCL.data];
%     Feature_name = [features.DrugCL.name];    
%     Feature_type = [features.DrugCL.type];    

        
%     factor = 2;
%     pos_labels = find(Syn_labels == 1);
%     neg_labels = find(Syn_labels == 0);
% 
%     for i = 1:20
%         fd_names = fopen(sprintf('Dream_Sample%d.names', i), 'w');
%         fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
%         for k = 1:numel(Feature_name)
%             fprintf(fd_names, '%s:\t%s.\n', Feature_name{k}, Feature_type{k});
%         end    
%         fclose(fd_names);
%         
%         neg_samples = randsample(neg_labels, numel(pos_labels)*factor);
%         samples = [pos_labels; neg_samples];
%         Sampled_Features = [num2cell(Syn_labels(samples, 1)), Feature_data(samples, :)];    
%         dlmcell(sprintf('Dream_Sample%d.data', i), Sampled_Features, 'delimiter', ',');    
%     end
        
        fd_names = fopen('Dream.names', 'w');
        fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
        for k = 1:numel(Feature_name)
            fprintf(fd_names, '%s:\t%s.\n', Feature_name{k}, Feature_type{k});
        end    
        fclose(fd_names);        
        dlmcell('Dream.data', [num2cell(Syn_labels), Feature_data], 'delimiter', ',');  

%%    

    empty = find(Syn_labels == 0);
    nonempty = find(Syn_labels == 1);
    rows = randsample(empty, round(1*numel(empty)));
    rows = union(rows, nonempty);
    
    Drug_Range = [12, 13, 22, 23, 26, 27, 28];
%     Drug_Range = 1:size(Feature_data, 2);
    
    Features = cell2mat(Feature_data(rows, Drug_Range));
%     Features = randn(numel(rows), numel(Drug_Range));
    Feature_names = Feature_name(Drug_Range);
    Labels = double(Syn_labels(rows));
    Label_names = cell(numel(Labels), 1);
    Label_names(:) = {'no'};
    Label_names(Labels == 1) = {'yes'};

    Feature_dataset = mat2dataset(Features, 'VarNames', Feature_names);

%%
    % Training set
    Xtrain = X(training(cv),:);
    Ytrain = Y(training(cv),:);
    Ytrain_val = Labels(training(cv));

    % Test set
    Xtest = X(test(cv),:);
    Ytest = Y(test(cv),:);
    Ytest_val = Labels(test(cv));

    disp('Training Set')
    tabulate(Ytrain)
    disp('Test Set')
    tabulate(Ytest)

%%
    DefaultTree = fitctree(Features, Labels, 'PredictorNames', Feature_names, 'MinLeaf', 50, 'PruneCriterion', 'error');
    view(DefaultTree,'Mode','Graph')
    
%%
    opts = statset('UseParallel',true);

    cost = [0 1
            1 0];
    tb = TreeBagger(50, Xtrain, Ytrain_val,'Method','Classification','OOBVarImp','On', 'cost', cost);

    [Y_tb, classifScore] = tb.predict(Xtest);    
    Y_tb = nominal(Y_tb);    

    % Compute the confusion matrix
    C_tb = confusionmat(Ytest_val, double(Y_tb)-1);
    % Examine the confusion matrix for each class as a percentage of the true class
    C_tb = bsxfun(@rdivide,C_tb,sum(C_tb,2)) * 100

     nnz(double(Y_tb) == 2)
    nnz(double(Y_tb) == 2 &  Ytest_val == 1) ./ nnz(Ytest_val)
    nnz(Ytest_val == 1)   
    
    [xx,yy,~,auc] = perfcurve(Ytest_val, classifScore(:,2), 1);
    figure;
    plot(xx,yy)
    xlabel('False positive rate');
    ylabel('True positive rate')
    title('ROC curve for ''yes'', predicted vs. actual response (Test Set)')
    text(0.5,0.25,{'TreeBagger with full feature set',strcat('Area Under Curve = ',num2str(auc))},'EdgeColor','k');

%     tb.OOBPermutedVarDeltaError  
    
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    h = bar(tb.OOBPermutedVarDeltaError);

    set(gca, 'XLim', [0, numel(tb.OOBPermutedVarDeltaError) + 1]);
    set(gca,'FontSize', 16, 'FontWeight','bold'); 
    xticklabel_rotate(1:numel(Feature_names), 45, Feature_names, 'FontSize', 16, 'FontWeight','bold', 'Interpreter', 'none');        
    set(h(1),'FaceColor', [0.5 0.5 0.5]);
%     print(fig, 'output/plots/feature_importance.eps', '-depsc2','-r300');
    
    
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    oobErrorBaggedEnsemble = oobError(tb);
    plot(oobErrorBaggedEnsemble, 'LineWidth', 2)
    set(gca,'FontSize', 14, 'FontWeight','bold'); 
    xlabel( 'Number of grown trees','FontSize', 16, 'FontWeight','bold' );
    ylabel('Out-of-bag classification error','FontSize', 16, 'FontWeight','bold' );   

    
%%
[~,idxvarimp] = sort(tb.OOBPermutedVarDeltaError, 'descend');
critfun = @(Xtr,Ytr,Xte,Yte) featureImp(Xtr,Ytr,Xte,Yte,'TreeBagger');
% The top 5 features determined in the previous step have been included,
% to reduce the number of combinations to be tried by sequentialfs
[fs,history] = sequentialfs(critfun, Xtrain, Ytrain,'keepin',idxvarimp(1:5),'Options',opts);
disp('Included features:');
disp(names(fs)');

%%
%     print(fig, 'output/plots/out-of-bag_error.eps', '-depsc2','-r300');
    
% %% Construct optimal Decision Tree and choose the right depth
%     leafs = logspace(1,2,10); % Generate minimum leaf occupancies for classification trees from 10 to 100, spaced exponentially apart
% 
%     rng('default')
%     N = numel(leafs);
%     err = zeros(N,1);
%     for n=1:N
%         t = fitctree(Features, Syn_labels, 'CrossVal', 'On', 'MinLeaf', leafs(n));
%         err(n) = kfoldLoss(t);
%     end
%     plot(leafs,err);
%     xlabel('Min Leaf Size');
%     ylabel('cross-validated error');
% 
% 
% %%
% DefaultTree = fitctree(Features, Syn_labels, 'PredictorNames', Feature_names, 'MinLeaf', 40);
% view(DefaultTree,'Mode','Graph')
% 
% OptimalTree = fitctree(Features, Syn_labels, 'minleaf', 20, 'PredictorNames', Feature_names);
% view(OptimalTree,'mode','graph')
% 
% resubOpt = resubLoss(OptimalTree);
% lossOpt = kfoldLoss(crossval(OptimalTree));
% resubDefault = resubLoss(DefaultTree);
% lossDefault = kfoldLoss(crossval(DefaultTree));
% resubOpt,resubDefault,lossOpt,lossDefault
% 
% %%
%     OptimalTree = fitctree(Features, Syn_labels, 'MinLeaf', 30, 'PredictorNames', Feature_names);
%     view(OptimalTree, 'mode','graph')
