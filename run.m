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
    % [ Mono ] = read_MonoTherapy(annotations, 'input/Dream/synergy/ch2_leaderBoard_monoTherapy.csv' );
    [ Mono, Pairs ] = read_MonoTherapy(annotations, 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv' );
    Pair_names = arrayfun(@(x) annotations.drugs.ChallengeName{x}, Pairs, 'UniformOutput', false);

    X = Mono.Drug_sensitivity';
    X(isnan(X) | isinf(X)) = 0;
%     [ii, jj, vv] = find(Mono.IC50);
%     sigma = std(vv);
%     vv = exp(vv / (2*sigma));
%     vv = vv ./ max(vv);
%     X = full(sparse(ii, jj, vv, size(annotations.cellLines, 1), size(annotations.drugs, 1)))';    
    X = normalize(X, 'dim', 2);
    
    for i = 1:size(Pairs, 1)
        for j = 1:size(annotations.cellLines, 1)
            Y(i, j) = (X(Pairs(i, 1), j) * X(Pairs(i, 2), j));
        end
    end
    Y = 1000* Y ./ (max(Y(:)));
    
    Z = Mono.Synergy;
    Z(Z < 0) = 0;
%     [ii, jj, vv] = find(Z);
%     ind = sub2ind(size(Z), ii, jj);
%     corr(vv, Y(ind))
    [cc, pval] = corr(Z(:), Y(:), 'type', 'Spearman')
%% TODO: Can we use MonoTherapy data as a source for computing D2D?
%     X = Mono.Drug_sensitivity; X(isnan(X)) = 0;
%     Y = double(logical(X));
%     XX = X'*X ./ (Y'*Y);
%     XX(isnan(XX)) = 0;
%     XX(XX > 0) = (nonzeros(XX) - min(nonzeros(XX))) / (max(nonzeros(XX)) - min(nonzeros(XX)));
%     SenSim = XX;
%     SenSim = SenSim - diag(diag(SenSim));
%     SenSim(SenSim < 0.9) = 0;
%     clustergram(SenSim, 'RowLabels', annotations.drugs.ChallengeName, 'ColumnLabels', annotations.drugs.ChallengeName, 'Linkage', 'average', 'ColorMap', colormap(flipud(redgreencmap())), 'OPTIMALLEAFORDER', true)


%% Drug-Drug Similarity Network
    D2D = Construct_D2D(annotations, ACSN);

    
%% Cellline-Celline Similarity Network
    % TODO: Which networks to use? Should we also use co-methyl and co-mut? How
    % to optimally combine using Mashup, GeneMANIA, or SNF?
    [C2C, NodeWeights] = Construct_C2C(annotations, ACSN, 'expression_only', true);


%% Compute topological gene signatures using cell type-specific interactomes
%     for beta = [0.01, 0.15, 0.5, 0.85, 0.99]
%         [ topological_gene_signature, topo_class_gene_idx ] = computeTopoSig( annotations, ACSN, NodeWeights, 'beta', beta);
%     end

    [ topological_gene_signature, topo_class_gene_idx ] = computeTopoSig( annotations, ACSN, NodeWeights, 'beta', 0.99);


%% Read LINCS dataset and compute transcriptional gene signatures
    [ transcriptional_gene_signature, transc_class_gene_idx ] = computeTransSig( annotations, ACSN );
    prop_Exp = propagate_Expression(D2D, C2C , transcriptional_gene_signature, 0.1);   
    % Cross-validate by masking 2 cell lines or 10 drugs

    
%% Expected signature of cancer drugs -- synergy among unctional classes
   [ ACSN_class_dir_weight, ACSN_class_corr] = computeTypicalTranSig( 'visualize', false  );
    
    
%% Compute Synergy scores
    % i) Integrate with Monotherapy to account for cell-type specific resistance    
    % ii) Synergy among functions: 

%     class_gene_idx = cellfun(@(genes) find(ismember(LINCS_genes, genes)), ACSN.class_genes, 'UniformOutput', false) ; % For transcriptional
    
    
   
    % Project signatures onto the functional space of ACSN classes
    fprintf('Projecting gene signatures onto the functional space of ACSN classes ... ');
    FunSig = cell(size(annotations.drugs, 1), size(annotations.cellLines, 1));
    for i = 1:size(annotations.drugs, 1)
        for j = 1:size(annotations.cellLines, 1)            
            FunSig{i, j} = cellfun(@(rows) topological_gene_signature{i, j}(rows), topo_class_gene_idx, 'UniformOutput', false);
        end
    end
    fprintf('done\n');
    
    
    Confidence_mat = nan(size(Pairs, 1), size(annotations.cellLines, 1));
    for pIdx = 1:size(Pairs, 1)
        fprintf('Pair %d/ %d\n', pIdx , size(Pairs, 1));
        d1 = Pairs(pIdx, 1);
        d2 = Pairs(pIdx, 2);
        for cIdx = 1:size(annotations.cellLines, 1)   
            fprintf('\t%d- Cell line %s\n', cIdx, annotations.cellLines.Sanger_Name{cIdx});
            S1 = FunSig{d1, cIdx};
            S2 = FunSig{d2, cIdx};
            if( isempty(S1) || isempty(S2) )
                continue;
            end
            empty_cells =  cellfun(@(x) isempty(x), S1) | cellfun(@(x) isempty(x), S2);
            S1(empty_cells) = [];
            S2(empty_cells) = [];

            FSyn = cell2mat(cellfun(@(s1, s2) (nanmean((abs(s1)+abs(s2))./(max(abs(s1), abs(s2)))) - 1)*(1 - abs(corr(s1, s2))), S1, S2, 'UniformOutput', false));
            synergy_score = nanmean(FSyn); % or sum(func_synergy_weight.*FSyn), where func_synergy_weight is the normalized directional importance of each class.

            % make sure how it affects scaling the synergy later
%             synergy_score = synergy_score * Mono.Drug_sensitivity(d1) * Mono.Drug_sensitivity(d2); 
            
            %TODO: Scaling by max is bad because it scales minor changes in
            %a class close to 1
            
            Confidence_mat(pIdx, cIdx) = synergy_score;
        end
    end

    vals = nonzeros(Confidence_mat);
    vals(isnan(vals) | isinf(vals)) = [];
    syn_threshold = quantile(vals, 0.8); % 50% of elements are predicted as synergic. This % is estimated from training data
    
    challenge_threshold = 30;    
    syn_factor = challenge_threshold/syn_threshold;
    Confidence_mat = Confidence_mat * syn_factor;
    

%% Export
    outPath = fullfile('output', 'predictions', 'leaderBoard');    
    exportResults( annotations, Pairs, Pair_names, Confidence_mat, outPath, 'synergy_threshold', 30);


