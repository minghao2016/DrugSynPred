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
%     [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' );
    [ Mono ] = read_singleMonoTherapy( annotations, 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv' );
    [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, 'input/Dream/synergy/ch1_train_combination_and_monoTherapy.csv' );

%%
clear ZZ
Z = Pair_synergy;
Z(isinf(Z)) = 0;
Z(Z < 0) = 0;
for i = 1:size(Z, 1)
    row = nonzeros(Z(i, :));
    val = median(row);
    if(~isempty(val))
        ZZ(i, 1) = val;    
    else
        ZZ(i, 1) = 0;
    end
%     ZZ(i) = max(Z(i, :));    
end    
ZZ(isnan(ZZ)) = 0;

%%
fd = fopen('/home/shahin/Dropbox/Dream/experiment/refs/GeneticInteractions/SynLethDB_human.txt', 'r');
C = textscan(fd, '%s %s');
fclose(fd);
SL_nodes = union(C{1}, C{2});
[~, ii] = ismember(C{1}, SL_nodes);
[~, jj] = ismember(C{2}, SL_nodes);
SL_net = sparse(ii, jj, 1, numel(SL_nodes), numel(SL_nodes));
SL_net = max(SL_net, SL_net');



% fd = fopen('/home/shahin/Dropbox/Dream/experiment/refs/GeneticInteractions/SL.sif', 'r');
% C = textscan(fd, '%s %s %s');
% SL_nodes = union(C{1}, C{3});
% [~, ii] = ismember(C{1}, SL_nodes);
% [~, jj] = ismember(C{3}, SL_nodes);
% SL_net = sparse(ii, jj, 1, numel(SL_nodes), numel(SL_nodes));
% SL_net = max(SL_net, SL_net');

% for i = 1:size(annotations.drugs, 1)
%     [~, SL_matches] = ismember(annotations.drugs.Target{i}, SL_nodes);
%     SL_matches(~SL_matches) = [];
%     SL_size(i) = numel(SL_matches);
% %     if(numel(SL_matches == 0))
% %         continue;
% %     end
%     SL_subnet = SL_net(SL_matches, :);
%     
%     genes{i} = {};
%     for j = 1:numel(SL_matches)
%         genes{i} = union(genes{i}, SL_nodes(find(SL_subnet(j, :))));
%     end    
% end

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
    N1 = union(t1_idx, N1);
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
    N2 = union(t2_idx, N2);
    N2_names = ACSN.vertex_names(N2);
    [~, SL_matches2] = ismember(N2_names, SL_nodes);
    SL_matches2(~SL_matches2) = [];

  
    Syn(i) = nnz(SL_net(SL_matches1, SL_matches2));
%     Syn(i) = max(numel(intersect(v1, SL_N2)), numel(intersect(v2, SL_N1)));
    
%     Syn(i) = numel(nonzeros(SL_net(SL_matches1, SL_matches2)));

    
%     
%     
%     Syn(i) = max(numel(intersect(v1, SL_N2)), numel(intersect(v2, SL_N1)));
    %     Syn(i) = SL_size(Pairs(i, 1)) + SL_size(Pairs(i, 2));    
end

x1 = tiedrank(ZZ);
x2 = tiedrank(Syn);
[c, p] = corr(x1, x2)









% Syn = zeros(size(Pairs, 1), 1);
% for i = 1:size(Pairs, 1)
%     t1 = annotations.drugs.Target{Pairs(i, 1)};
%     [~, t1_idx] = ismember(t1, ACSN.vertex_names); 
%     t1_idx(~t1_idx) = [];
%     t1(~t1_idx) = [];
%     t1_subnet = ACSN.A(t1_idx, :);
%     N1 = [];
%     for j = 1:numel(t1_idx)
%         N1 = union(N1, find(t1_subnet(j, :)));
%     end
%     v1 = union(t1_idx, N1);
%     v1_names = ACSN.vertex_names(v1);
%     [~, SL_matches1] = ismember(v1_names, SL_nodes);
%     SL_matches1(~SL_matches1) = [];
%     
%     
%     SL_subnet = SL_net(SL_matches, :);    
%     SL_N1_genes = {};
%     for j = 1:numel(SL_matches)
%         SL_N1_genes = union(SL_N1_genes, SL_nodes(find(SL_subnet(j, :))));
%     end 
%     [~, SL_N1] = ismember(SL_N1_genes, ACSN.vertex_genes);
%     SL_N1(~SL_N1) = [];
%     
%     
%     t2 = annotations.drugs.Target{Pairs(i, 2)};
%     [~, t2_idx] = ismember(t2, ACSN.vertex_names); 
%     t2_idx(~t2_idx) = [];
%     t2(~t2_idx) = [];
%     t2_subnet = ACSN.A(t2_idx, :);
%     N2 = [];
%     for j = 2:numel(t2_idx)
%         N2 = union(N2, find(t2_subnet(j, :)));
%     end
%     v2 = union(t2_idx, N2);
%     v2_names = ACSN.vertex_names(v2);
%     [~, SL_matches2] = ismember(v2_names, SL_nodes);
%     SL_matches2(~SL_matches2) = [];
% 
%     SL_subnet = SL_net(SL_matches, :);
%     SL_N2_genes = {};
%     for j = 2:numel(SL_matches)
%         SL_N2_genes = union(SL_N2_genes, SL_nodes(find(SL_subnet(j, :))));
%     end 
%     [~, SL_N2] = ismember(SL_N2_genes, ACSN.vertex_genes);
%     SL_N2(~SL_N2) = [];
%     
%     Syn(i) = max(numel(intersect(v1, SL_N2)), numel(intersect(v2, SL_N1)));
%     
% %     Syn(i) = numel(nonzeros(SL_net(SL_matches1, SL_matches2)));
% 
%     
% %     
% %     
% %     Syn(i) = max(numel(intersect(v1, SL_N2)), numel(intersect(v2, SL_N1)));
%     %     Syn(i) = SL_size(Pairs(i, 1)) + SL_size(Pairs(i, 2));    
% end
% 





% Syn = zeros(size(Pairs, 1), 1);
% for i = 1:size(Pairs, 1)
%     t1 = annotations.drugs.Target{Pairs(i, 1)};
%     [~, t1_idx] = ismember(t1, ACSN.vertex_names); 
%     t1_idx(~t1_idx) = [];
%     t1(~t1_idx) = [];
%     t1_subnet = ACSN.A(t1_idx, :);
%     N1 = [];
%     for j = 1:numel(t1_idx)
%         N1 = union(N1, find(t1_subnet(j, :)));
%     end
%     v1 = union(t1_idx, N1);
%     v1_names = ACSN.vertex_names(v1);
%     [~, SL_matches1] = ismember(v1_names, SL_nodes);
%     SL_matches1(~SL_matches1) = [];
%     [~, v1_SL_idx] = ismember(t1, SL_nodes);
%     v1_SL_idx(~v1_SL_idx) = [];
%     if(isempty(v1_SL_idx))
%         continue;
%     end
%     
%     
%     
%     t2 = annotations.drugs.Target{Pairs(i, 2)};
%     [~, t2_idx] = ismember(t2, ACSN.vertex_names); 
%     t2_idx(~t2_idx) = [];
%     t2(~t2_idx) = [];
%     t2_subnet = ACSN.A(t2_idx, :);
%     N2 = [];
%     for j = 2:numel(t2_idx)
%         N2 = union(N2, find(t2_subnet(j, :)));
%     end
%     v2 = union(t2_idx, N2);
%     v2_names = ACSN.vertex_names(v2);
%     [~, SL_matches2] = ismember(v2_names, SL_nodes);
%     SL_matches2(~SL_matches2) = [];
%     [~, v2_SL_idx] = ismember(t2, SL_nodes);
%     v2_SL_idx(~v2_SL_idx) = [];
%     if(isempty(v2_SL_idx))
%         continue;
%     end
%     
%     Syn(i) = max(numel(nonzeros(SL_net(v1_SL_idx, SL_matches2))), numel(nonzeros(SL_net(v2_SL_idx, SL_matches1))));
% end







%%
alpha = 0.85;
n = size(ACSN.A, 1);
e = ones(n, 1) ./ n;
lambda = max(eigs(ACSN.A));
A_org = ACSN.A ./ lambda;

frows = [];
for i = 1:numel(ACSN.class_names)
    [~, crows] = ismember(ACSN.class_genes{i}, ACSN.vertex_genes);
    frows = union(frows, crows);
end
frows(~frows) = [];

for i = 86:size(Pairs, 1)
    fprintf('%d- %s vs %s \n', i, annotations.drugs.ChallengeName{Pairs(i, 1)}, annotations.drugs.ChallengeName{Pairs(i, 2)});
    t1 = annotations.drugs.Target{Pairs(i, 1)};
    t2 = annotations.drugs.Target{Pairs(i, 2)};
    
    [~, t1_rows] = ismember(t1, ACSN.vertex_genes); t1_rows(t1_rows == 0) = [];
    [~, t2_rows] = ismember(t2, ACSN.vertex_genes); t2_rows(t2_rows == 0) = [];
   
    A = A_org;
    A(t1_rows, :) = 0;
    A(:, t1_rows) = 0;
    sig{i, 1} = (1-alpha).*(eye(n) - alpha*A)\e;
    
    A = A_org;
    A(t2_rows, :) = 0;
    A(:, t2_rows) = 0;
    sig{i, 2} = (1-alpha).*(eye(n) - alpha*A)\e;
   
    A(t1_rows, :) = 0;
    A(:, t1_rows) = 0;
    sig{i, 3} = (1-alpha).*(eye(n) - alpha*A)\e;
    
    Delta = log(geomean([sig{i, 1}, sig{i, 2}], 2) ./ sig{i, 3});
    fDelta = Delta(frows);    
    
    M(i) = nansum(Delta);
    
    
    M2(i) = nansum(fDelta);
    
end


for i = 1:size(Pairs, 1)
    fprintf('%d- %s vs %s \n', i, annotations.drugs.ChallengeName{Pairs(i, 1)}, annotations.drugs.ChallengeName{Pairs(i, 2)});

    Delta = log(max([sig{i, 1}, sig{i, 2}], [], 2) ./ sig{i, 3});
%     Delta = sig{i, 3};

    fDelta = Delta(frows);    

    Delta(isnan(Delta) |isinf(Delta)) = [];
    fDelta(isnan(fDelta) |isinf(fDelta)) = [];
    
    M(i) = nansum(Delta);
    M2(i) = nansum(fDelta);   
end



[cc, pval] = corr(M2, ZZ)
%%

    X = Mono.Drug_sensitivity;
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
    
    Z = Pair_synergy;
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


