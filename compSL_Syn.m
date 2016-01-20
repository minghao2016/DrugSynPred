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
    [ VertSyn_Topo ] = computeVertSyn_topo( annotations, interactome, Pairs );

%%
    [ VertSyn_Expr, VertSyn_Expr_Rest ] = computeVertSyn_expr( annotations, Pairs );
    
    %%
    Drug_Effect = cell2mat(Mono.EMax);
    Drug_Effect(isnan(Drug_Effect)) = 100;
    Drug_Effect = 100 - Drug_Effect; % effectiveness of drug in cell-lines
    
    [ Horiz_compSL_Syn ] = computeHorSyn_DrugEffect( annotations, Pairs, Drug_Effect );    
    
%%
    % Only works AFTER computing Synergy scores. If we move it to the
    % beginning results are horrible


%     D2D = D2D_combined(annotations, interactome); % Construct_D2D(annotations, interactome);
%     [C2C, NodeWeights] = Construct_C2C(annotations, interactome, 'expression_only', true);    
%     X = cell2mat(fillinBlanks(Mono.EMax, C2C, D2D));
    X = cell2mat(Mono.EMax);
    
    X(isempty(X) | isnan(X) | isinf(X)) = 0; 
    X = normalize(X, 'dim', 1);    

    TT = Horiz_compSL_Syn;
    TT(TT > 1e-1) = 0;
    V{1} = full(spfun(@(x) -log10(x), TT));

%     V{1} = -log10(Horiz_compSL_Syn);
    
    Z = zeros(size(Pairs, 1), size(annotations.cellLines, 1));
    for i = 1:size(Pairs, 1)
        for j = 1:size(annotations.cellLines, 1)
            Z(i, j) = V{1}(i) * max([X(Pairs(i, 1), j), X(Pairs(i, 2), j)], [], 2);
        end
    end    
%     Z = Modified_zscore(Z')';
%     Z(isnan(Z)) = 0;
%     Z = Z ./ max(nonzeros(Z));

%     [z, z_idx] = sort(nonzeros(Z));
%     z_cut = cut(z, 'selection_method', 'participation_ratio_stringent');
%     z_threshold = z(z_cut);

    z_threshold = 30;
    Z = 30*Z ./ prctile(nonzeros(Z), 90);
    ZZ = Z > z_threshold;


    outPath = fullfile('output', 'predictions', 'leaderBoard');        
    exportResults( annotations, Pairs, Pair_names, Z, outPath, 'synergy_threshold', z_threshold);

    

    