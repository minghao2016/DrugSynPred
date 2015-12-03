function [Pairs, Pair_names, Pair_synergy, Pair_quality] = readPairs( annotations, fname)

    fprintf('Importing drug pairs info ...\n');
        fprintf('\t%s\n', fname);

    T = readtable(fname, 'TreatAsEmpty', 'NA', 'Format', '%s%s%s%f%f%f%f%f%f%f%f%f%f%s', 'Delimiter', ',');
    [~, Drug_idx_A] = ismember(T.COMPOUND_A, annotations.drugs.ChallengeName);
    [~, Drug_idx_B] = ismember(T.COMPOUND_B, annotations.drugs.ChallengeName);
    [~, CL_idx] = ismember(T.CELL_LINE, annotations.cellLines.Sanger_Name);

    Pairs = unique([Drug_idx_A, Drug_idx_B], 'rows');
    Pair_names = arrayfun(@(x) annotations.drugs.ChallengeName{x}, Pairs, 'UniformOutput', false);


    Pair_synergy = -inf(size(Pairs, 1), size(annotations.cellLines, 1));
    Pair_quality = -inf(size(Pairs, 1), size(annotations.cellLines, 1));

    for i = 1:size(T, 1)
        if(~isnan(T.SYNERGY_SCORE(i)))
            pair_mask = (Pairs(:, 1) == Drug_idx_A(i) & Pairs(:, 2) == Drug_idx_B(i));
            Pair_synergy(pair_mask, CL_idx(i)) = T.SYNERGY_SCORE(i);
        end
        if(~isnan(T.QA(i)))
            Pair_quality(Drug_idx_A(i), Drug_idx_B(i)) = T.QA(i);
        end
    end        
      
end

