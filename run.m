clear
addpath(genpath('code'));
warning('off','all');    

annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);
sorted_CL = sort(annotations.cellLines.Sanger_Name);

[~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
annotations.cellLines = annotations.cellLines(CL_perm, :);

experiment_type = 'leaderBoard';

%% Load signaling interactome from ACSN together with 55 functional classes relevant to cancer
% TODO: To expand nr not to expand?
% TODO\: Edge directions
[ ACSN ] = import_ACSN();


%% Read Monotherapy data and impute missing values based on Dual Layer method
% [ Mono ] = read_MonoTherapy(annotations, 'input/Dream/synergy/ch2_leaderBoard_monoTherapy.csv' );
[ Mono, Pairs ] = read_MonoTherapy(annotations, 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv' );
Pair_names = arrayfun(@(x) annotations.drugs.ChallengeName{x}, Pairs, 'UniformOutput', false);

%% TODO: Can we use MonoTherapy data as a source for computing D2D?
% X = Mono.Drug_sensitivity; X(isnan(X)) = 0;
% Y = double(logical(X));
% XX = X'*X ./ (Y'*Y);
% XX(isnan(XX)) = 0;
% XX(XX > 0) = (nonzeros(XX) - min(nonzeros(XX))) / (max(nonzeros(XX)) - min(nonzeros(XX)));
% SenSim = XX;
% SenSim = SenSim - diag(diag(SenSim));
% SenSim(SenSim < 0.9) = 0;
% clustergram(SenSim, 'RowLabels', annotations.drugs.ChallengeName, 'ColumnLabels', annotations.drugs.ChallengeName, 'Linkage', 'average', 'ColorMap', colormap(flipud(redgreencmap())), 'OPTIMALLEAFORDER', true)


%% Drug-Drug Similarity Network
% TODO: Gene targets for DNA, Methylation, ... targeting drugs
targetD2D = Construct_D2D(ACSN, annotations);
stitchD2D = Construct_stitchD2D();

D2D = targetD2D + stitchD2D ; 

%% Cellline-Celline Similarity Network
% TODO: Which networks to use? Should we also use co-methyl and co-mut? How
% to optimally combine using Mashup, GeneMANIA, or SNF?
C2C = Construct_C2C(annotations, 'expression_only', true);

% 
% %% Read Leadership board
%     Leadership = readtable('input/Dream/submission/ch2/leadership/synergy_matrix.csv', 'ReadRowNames', true);
%     pair_names = cellfun(@(x) strsplit(x, '.'), Leadership.Properties.RowNames, 'UniformOutput', false);
%     pair_idx = cell(numel(pair_names), 1);
%     for i = 1:numel(pair_names)
%         pair_idx{i}{1} = find(strcmp(pair_names{i}{1}, annotations.drugs.ChallengeName));
%         pair_idx{i}{2} = find(strcmp(pair_names{i}{2}, annotations.drugs.ChallengeName));        
%     end
%     sorted_CL = Leadership.Properties.VariableNames;

%% Read LINCS dataset
    % TODO: Draft file! CHECK CHECK CHECK, to make sure we selected
    % the best drugs to assay
    LINCS_ds = parse_gct('input/LINCS/final/LINCS_subset.gct');
   
    LINCS_celllines = LINCS_ds.cdesc(:, 1);
    LINCS_celllines(strcmp(LINCS_celllines, 'BT20')) = {'BT-20'};
    LINCS_celllines(strcmp(LINCS_celllines, 'HT29')) = {'HT-29'};    
    LINCS_celllines(strcmp(LINCS_celllines, 'MDAMB231')) = {'MDA-MB-231'};    
    
    LINCS_drugs = LINCS_ds.cdesc(:, 7);
    LINCS_expression_matrix = LINCS_ds.mat; % TODO: Should we column normalize to ensure constant transcriptional activity for each drug?
    LINCS_expression_within_groups = zeros(numel(ACSN.class_names), size(LINCS_expression_matrix, 2));
    
    
    for g = 1:numel(ACSN.class_genes)
        [~, rows] = ismember(ACSN.class_genes{g}, LINCS_ds.rdesc(:, 7));
        rows(rows == 0) = [];
        LINCS_expression_within_groups(g, :) = arrayfun(@(col) mean(LINCS_expression_matrix(rows, col)), 1:size(LINCS_expression_matrix, 2));        
    end
    
    
    
    [~, cl_idx] = ismember(LINCS_celllines, annotations.cellLines.Sanger_Name);
    % TODO: Check to make sure all ID mappings are correct
    Dream2LINCS= readtable('./input/LINCS/final/preliminary_mapping.csv');
    
    Expr_DS = cell(size(annotations.drugs, 1), size(annotations.cellLines, 1));
    for i = 1:size(LINCS_expression_matrix, 2)        
        rows = find(ismember(Dream2LINCS.ID, LINCS_drugs{i}));        
        % TODO: Should we use aggregated scores in groups (probably), or all
        % genes without grouping (unlikely)?
        Expr_DS(rows, cl_idx(i)) = {LINCS_expression_matrix(:, i)};
        % *** OR ***
%         Expr_DS(rows, cl_idx(i)) = {LINCS_expression_within_groups(:, i)};
    end
    

    %Propagate the expression per cell line and per drug to impute missing
    %values
    alpha = 0.1;
    [D2D_prop_Exp,C2C_prop_Exp] = propagate_Expression(D2D, C2C , Expr_DS, alpha);
    
    % TODO: Use 2 Layer method to impute missing values

    
%% Compute Synergy scores
    % TODO: Synergy prediction
    % 1- Similarity
    % a- correlation between all genes in the pair of expression signatures
    % b- correlation between pathway genes in the pair of expression signatures
    % c- correlation between average of genes in each pathway in the pair of expression signatures
 
    % 2- Dissimilarity
    % a- Manhattan distance between all genes in the pair of expression signatures
    % b- Manhattan distance between pathway genes in the pair of expression signatures
    % c- Manhattan distance between average of genes in each pathway in the pair of expression signatures
    
    % 3- Drug topological signatures using TS networks
    
    % 4- Integrate with Monotherapy to account for cell-type specific
    % resistance
    
    % 5- Synergy among functions: 
    % a) nondirectionl (RWR over ACSN)
    
    Confidence_mat = nan(size(Pairs, 1), size(annotations.cellLines, 1));
    for pIdx = 1:size(Pairs, 1)
        d1 = Pairs(pIdx, 1);
        d2 = Pairs(pIdx, 2);
        for cIdx = 1:size(annotations.cellLines, 1)
            Confidence_mat(pIdx, cIdx) = 2 * rand(1);
        end
    end
    

%% Export
    synergy_threshold = 1; % TODO: What is the optimal threshold??
    
    fd_syn = fopen(fullfile('output', 'predictions', experiment_type, 'synergy_matrix.csv'), 'w');
    fd_conf = fopen(fullfile('output', 'predictions', experiment_type, 'confidence_matrix.csv'), 'w');
    for cIdx = 1:size(annotations.cellLines, 1)
        fprintf(fd_syn, ',%s', sorted_CL{cIdx});
        fprintf(fd_conf, ',%s', sorted_CL{cIdx});
    end
    fprintf(fd_syn, '\n');
    fprintf(fd_conf, '\n');    
    
    for pIdx = 1:size(Pairs, 1)
        fprintf(fd_syn, '%s.%s', Pair_names{pIdx, 1}, Pair_names{pIdx, 2});
        fprintf(fd_conf, '%s.%s', Pair_names{pIdx, 1}, Pair_names{pIdx, 2});
        for cIdx = 1:size(annotations.cellLines, 1)
            fprintf(fd_syn, ',%d', Confidence_mat(pIdx, cIdx) > synergy_threshold);
            fprintf(fd_conf, ',%f', Confidence_mat(pIdx, cIdx));
        end
        if(pIdx ~= size(Pairs, 1))
            fprintf(fd_syn, '\n');
            fprintf(fd_conf, '\n');    
        end
    end    
    fclose(fd_syn);
    fclose(fd_conf);


