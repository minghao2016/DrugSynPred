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
    [ Mono ] = read_allMonoTherapy( annotations, 'input/Dream/synergy/' );

%     [ Mono ] = read_singleMonoTherapy( annotations, fname );
    [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname );
