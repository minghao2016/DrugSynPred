    clear
    addpath(genpath('Classification'));
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
    
%%    
    Mono.EMax = fillinBlanks(Mono.EMax, C2C, D2D);
    Mono.H = fillinBlanks(Mono.H, C2C, D2D);
    Mono.IC50 = fillinBlanks(Mono.IC50, C2C, D2D);
    Mono.Drug_sensitivity = fillinBlanks(Mono.Drug_sensitivity, C2C, D2D);

     

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
    MW_data(:, 1) = annotations.drugs.MW(di);
    MW_data(:, 2) = annotations.drugs.MW(dj);
    MW_name = {'MW_1', 'MW_2'};

%     HBA_data(:, 1) = annotations.drugs.HBA(di);
%     HBA_data(:, 2) = annotations.drugs.HBA(dj);
%     HBA_name = {'HBA_1', 'HBA_2'};
% 
%     HBD_data(:, 1) = annotations.drugs.HBD(di);
%     HBD_data(:, 2) = annotations.drugs.HBD(dj);
%     HBD_name = {'HBD_1', 'HBD_2'};

    cLogP_data(:, 1) = annotations.drugs.cLogP(di);
    cLogP_data(:, 2) = annotations.drugs.cLogP(dj);
    cLogP_name = {'cLogP_1', 'cLogP_2'};

    Lipinski_data(:, 1) = annotations.drugs.Lipinski(di);
    Lipinski_data(:, 2) = annotations.drugs.Lipinski(dj);
    Lipinski_name = {'Lipinski_1', 'Lipinski_2'};

%     features.Drug.data = num2cell([MW_data, HBA_data, HBD_data, cLogP_data, Lipinski_data]);
%     features.Drug.name = [MW_name, HBA_name, HBD_name, cLogP_name, Lipinski_name]';
%     features.Drug.type = repmat({'continuous'}, numel(features.Drug.name), 1);

    features.Drug.data = num2cell([MW_data, cLogP_data, Lipinski_data]);
    features.Drug.name = [MW_name, cLogP_name, Lipinski_name]';
    features.Drug.type = repmat({'continuous'}, numel(features.Drug.name), 1);
%%  Drug-CellLine features

% 	IC50_data(:, 1) = [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]';
% 	IC50_data(:, 2) = [Mono.IC50{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
%     IC50_data = -log10(IC50_data);
%     IC50_data = (nanmax(IC50_data(:)) - IC50_data) ./ (nanmax(IC50_data(:)) - nanmin(IC50_data(:)));
% 	IC50_name = {'IC50_1', 'IC50_2'};
% 
% 	H_data(:, 1) = [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]';
% 	H_data(:, 2) = [Mono.H{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
%     H_data = (nanmax(H_data(:)) - H_data) ./ (nanmax(H_data(:)) - nanmin(H_data(:)));    
% 	H_name = {'H_1', 'H_2'};
    
	EMax_data(:, 1) = [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]';
	EMax_data(:, 2) = [Mono.EMax{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
    EMax_data = 1 - (EMax_data ./ 100);
	EMax_name = {'EMax_1', 'EMax_2'};    

	Drug_sensitivity_data(:, 1) = [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], di, ck)}]';
	Drug_sensitivity_data(:, 2) = [Mono.Drug_sensitivity{sub2ind([size(annotations.drugs, 1), size(annotations.cellLines, 1)], dj, ck)}]';
    Drug_sensitivity_data = (nanmax(Drug_sensitivity_data(:)) - Drug_sensitivity_data) ./ (nanmax(Drug_sensitivity_data(:)) - nanmin(Drug_sensitivity_data(:)));    
	Drug_sensitivity_name = {'Drug_sensitivity_1', 'Drug_sensitivity_2'};
    
%     features.DrugCL.data = num2cell([IC50_data, H_data, EMax_data, Drug_sensitivity_data]);
%     features.DrugCL.name = [IC50_name, H_name, EMax_name, Drug_sensitivity_name]';
%     features.DrugCL.type = repmat({'continuous'}, numel(features.DrugCL.name), 1);

    
%     [CC, pval] = corr(cell2mat(features.DrugCL.data), Syn_labels)
% 
%     [CC, pval] = corr(cell2mat(features.DrugCL.data), Syn_vv)


    features.DrugCL.data = num2cell([EMax_data, Drug_sensitivity_data]);
    features.DrugCL.name = [EMax_name, Drug_sensitivity_name]';
    features.DrugCL.type = repmat({'continuous'}, numel(features.DrugCL.name), 1);

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
    
    % Combine
    features.DrugDrug.data = [features.DrugDrug.data, num2cell(Horiz_compSL_Syn(Syn_ii))];
    features.DrugDrug.name = [features.DrugDrug.name; 'HorizSyn_Effect'];
    features.DrugDrug.type = [features.DrugDrug.type; 'continuous'];    
    
    
%% Drug-Drug-CL
%     % a) SL_sim only using expressed genes
%     
%     % b) Similarity of gene expression: 
%     % 1) Read LINC dataset
%     % 2) Compute Sim/Dissimilarity (different measures)
%     % 3) Project onto Census genes or ACSN functional classes before sim/dissim compitations
     
%% Combine all features
    Drug_perm = zeros(size(features.Drug.data, 2), 1);
    Drug_perm(1:2:numel(Drug_perm)) = 2:2:numel(Drug_perm);
    Drug_perm(2:2:numel(Drug_perm)) = 1:2:numel(Drug_perm);

    DrugCL_perm = zeros(size(features.DrugCL.data, 2), 1);
    DrugCL_perm(1:2:numel(DrugCL_perm)) = 2:2:numel(DrugCL_perm);
    DrugCL_perm(2:2:numel(DrugCL_perm)) = 1:2:numel(DrugCL_perm);
    
    Feature_data = [features.Drug.data, features.DrugCL.data, features.DrugDrug.data];
    Feature_data = [Feature_data; features.Drug.data(:, Drug_perm), features.DrugCL.data(:, DrugCL_perm), features.DrugDrug.data];

    
    Feature_name = [features.Drug.name; features.DrugCL.name; features.DrugDrug.name];    
    Feature_type = [features.Drug.type; features.DrugCL.type; features.DrugDrug.type]; 
    Feature_labels = double([Syn_labels;Syn_labels]);

    
    row_perm = randperm(size(Feature_data, 1), size(Feature_data, 1));
    Feature_data = Feature_data(row_perm, :);
    Feature_labels = Feature_labels(row_perm);
    
    Feature_Label_names = cell(numel(Feature_labels), 1);
    Feature_Label_names(:) = {'no'};
    Feature_Label_names(Feature_labels == 1) = {'yes'};

    
    
%%

%     Feature_data = [features.DrugCL.data, features.DrugDrug.data];
%     Feature_name = [features.DrugCL.name; features.DrugDrug.name];    
%     Feature_type = [features.DrugCL.type; features.DrugDrug.type];  

    
%     Feature_data = [features.CL.data, features.Drug.data, features.DrugCL.data, features.DrugDrug.data];
%     Feature_name = [features.CL.name; features.Drug.name; features.DrugCL.name; features.DrugDrug.name];    
%     Feature_type = [features.CL.type; features.Drug.type; features.DrugCL.type; features.DrugDrug.type];    
%     Feature_data = [features.DrugCL.data];
%     Feature_name = [features.DrugCL.name];    
%     Feature_type = [features.DrugCL.type];    


%%
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
        
%         fd_names = fopen('Dream.names', 'w');
%         fprintf(fd_names, 'Synergy.\n\nSynergy:\t0,1.\n');
%         for k = 1:numel(Feature_name)
%             fprintf(fd_names, '%s:\t%s.\n', Feature_name{k}, Feature_type{k});
%         end    
%         fclose(fd_names);        
%         dlmcell('Dream.data', [num2cell(Feature_labels), Feature_data], 'delimiter', ',');  

%%  Sub-sample

    empty = find(Feature_labels == 0);
    nonempty = find(Feature_labels == 1);
    rows = randsample(empty, round(1*numel(empty)));
    rows = union(rows, nonempty);
    
%     col_Range = [12, 13, 22, 23, 26, 27, 28];
    col_Range = 1:size(Feature_data, 2);
    
    X = cell2mat(Feature_data(rows, col_Range));
    Y = Feature_Label_names(rows);



%%
    cv = cvpartition(length(Y),'holdout',0.40);

    % Training set
    Xtrain = X(training(cv),:);
    Ytrain = Y(training(cv),:);
    Ytrain_val = Feature_labels(training(cv));

    % Test set
    Xtest = X(test(cv),:);
    Ytest = Y(test(cv),:);
    Ytest_val = Feature_labels(test(cv));

    disp('Training Set')
    tabulate(Ytrain)
    disp('Test Set')
    tabulate(Ytest)

    cost = [0 1
            1 0];    
        
%% SVM
    opts = statset('MaxIter',30000);
    % Train the classifier
    svmStruct = svmtrain(Xtrain,Ytrain,'kernel_function','rbf','kktviolationlevel',0.1,'options',opts);

    % Make a prediction for the test set
    [Y_svm, classifScore_svm]  = svmclassify(svmStruct,Xtest);
    C_svm = confusionmat(Ytest,Y_svm);
    % Examine the confusion matrix for each class as a percentage of the true class
    C_svm = bsxfun(@rdivide,C_svm,sum(C_svm,2)) * 100
    sum(diag(C_svm))

    [xx,yy,~,auc] = perfcurve(Ytest_val, -classifScore_svm, 1);
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    plot(xx,yy)
    xlabel('False positive rate');
    ylabel('True positive rate')
    title('ROC curve for ''yes'', predicted vs. actual response (Test Set)')
    text(0.5,0.25,{'TreeBagger with full feature set',strcat('Area Under Curve = ',num2str(auc))},'EdgeColor','k');
    print(fig, 'output/plots/SVM_AUC.eps', '-depsc2','-r300');    

%%
%     DefaultTree = fitctree(X, Y, 'PredictorNames', Feature_name, 'MinLeaf', 100, 'PruneCriterion', 'error', 'Cost', cost);
%     view(DefaultTree,'Mode','Graph')
    
%%
    opts = statset('UseParallel',true);


    tb = TreeBagger(500, Xtrain, Ytrain,'Method','Classification','OOBVarImp','On', 'cost', cost);

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
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    plot(xx,yy)
    xlabel('False positive rate');
    ylabel('True positive rate')
    title('ROC curve for ''yes'', predicted vs. actual response (Test Set)')
    text(0.5,0.25,{'TreeBagger with full feature set',strcat('Area Under Curve = ',num2str(auc))},'EdgeColor','k');
    print(fig, 'output/plots/RF_AUC.eps', '-depsc2','-r300');

    
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    h = bar(tb.OOBPermutedVarDeltaError);

    set(gca, 'XLim', [0, numel(tb.OOBPermutedVarDeltaError) + 1]);
    set(gca,'FontSize', 16, 'FontWeight','bold'); 
    xticklabel_rotate(1:numel(Feature_name), 45, Feature_name, 'FontSize', 16, 'FontWeight','bold', 'Interpreter', 'none');        
    set(h(1),'FaceColor', [0.5 0.5 0.5]);
    print(fig, 'output/plots/RF_feature_importance.eps', '-depsc2','-r300');
    
    fig = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1], 'PaperPositionMode', 'auto');  
    oobErrorBaggedEnsemble = oobError(tb);
    plot(oobErrorBaggedEnsemble, 'LineWidth', 2)
    set(gca,'FontSize', 14, 'FontWeight','bold'); 
    xlabel( 'Number of grown trees','FontSize', 16, 'FontWeight','bold' );
    ylabel('Out-of-bag classification error','FontSize', 16, 'FontWeight','bold' );
    print(fig, 'output/plots/RF_OOBError.eps', '-depsc2','-r300');
    

% %%
% [~,idxvarimp] = sort(tb.OOBPermutedVarDeltaError, 'descend');
% critfun = @(Xtr,Ytr,Xte,Yte) featureImp(Xtr,Ytr,Xte,Yte,'TreeBagger');
% % The top 5 features determined in the previous step have been included,
% % to reduce the number of combinations to be tried by sequentialfs
% [fs,history] = sequentialfs(critfun, Xtrain, Ytrain_val,'keepin',idxvarimp(1:5),'Options',opts);
% disp('Included features:');
% disp(Feature_name(fs)');

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
