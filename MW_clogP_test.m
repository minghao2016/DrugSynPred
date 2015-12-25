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
%     [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' ); %TODO: How to optimally combine replicates?

    [ Mono ] = read_singleMonoTherapy( annotations, fname );
    [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );

%%    
    X = Pair_synergy;
    X(isinf(X)) = 0;
    for i = 1:size(annotations.drugs, 1)
        pair_rows = [find(Pairs(:, 1) == i); find(Pairs(:, 2) == i)];
        drug_syns{i} = [];
        for j = 1:numel(pair_rows)
            drug_syns{i} = [drug_syns{i}; nonzeros(X(pair_rows(j), :))];
        end
    end

    ii = find(~cellfun(@(x) isempty(x), drug_syns));
    stats = cell2mat(cellfun(@(x) [nnz(x>20), nnz(x > 30), nnz(x >40), min(x), median(x), mean(x), max(x)], drug_syns(ii), 'UniformOutput', false)');
    
    MW = annotations.drugs.MW(ii);
    clogP = annotations.drugs.cLogP(ii);
    Lipinski = annotations.drugs.Lipinski(ii);
    

    
    Results = [stats, MW, clogP, Lipinski];
    Results(isnan(Results(:, 5)), :) = [];
    Results(:, 1:6) = Modified_zscore(Results(:, 1:6));
    
    Results(isnan(Results)) = 100*rand(1);
    [cc, pp] = corr(Results)
    
%     clear Results_table
%     Results_table(1, :) = {'Min_Syn', 'Median_Syn', 'Avg_Syn', 'Max_Syn', 'MW', 'cLogP', 'Lipinski'};
%     Results_table(2:size(Results, 1)+1, 1:size(Results, 2)) = num2cell(Results);
%     dlmcell('Drug_Promiscuity.tsv', Results_table);
    
%     for i = 1:4
%         for j = 5:7
%             rows = find(~isnan(Results(:, j)));
%             [cc, pp] = corr(Results(rows, i), Results(rows, j))
%             [cc, pp] = corr(
%             fprintf('%s vs %s: cc = %f, pp = %f\n', Results_table(1, i), Results_table(1, j)
%         end
%     end
%     [CC, PP] = corr(Results);