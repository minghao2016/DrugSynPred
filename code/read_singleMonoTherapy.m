function [ Mono ] = read_singleMonoTherapy( annotations, fname )


    fprintf('Importing Mono therapy data ...\n');
    % Read Mono therapy
    Mono.IC50 = inf(size(annotations.cellLines, 1), size(annotations.drugs, 1));
    Mono.EMax = inf(size(annotations.cellLines, 1), size(annotations.drugs, 1));
    Mono.H = nan(size(annotations.cellLines, 1), size(annotations.drugs, 1));
    Mono.Max_C = nan(size(annotations.cellLines, 1), size(annotations.drugs, 1));

    fprintf('\t%s\n', fname);
    T = readtable(fname, 'TreatAsEmpty', 'NA', 'Format', '%s%s%s%f%f%f%f%f%f%f%f%f%f%s', 'Delimiter', ',');
    [~, CL_idx] = ismember(T.CELL_LINE, annotations.cellLines.Sanger_Name);
    [~, Drug_idx_A] = ismember(T.COMPOUND_A, annotations.drugs.ChallengeName);
    [~, Drug_idx_B] = ismember(T.COMPOUND_B, annotations.drugs.ChallengeName);
    for i = 1:size(T, 1)

        if( T.IC50_A(i) < Mono.IC50(CL_idx(i), Drug_idx_A(i)) && Mono.EMax(CL_idx(i), Drug_idx_A(i)) > T.Einf_A(i))
            Mono.IC50(CL_idx(i), Drug_idx_A(i)) = T.IC50_A(i);
            Mono.EMax(CL_idx(i), Drug_idx_A(i)) = T.Einf_A(i);
            Mono.H(CL_idx(i), Drug_idx_A(i)) = T.H_A(i);        
            Mono.Max_C(CL_idx(i), Drug_idx_A(i)) = T.MAX_CONC_A(i);        
        end
        if( T.IC50_B(i) < Mono.IC50(CL_idx(i), Drug_idx_B(i)) && Mono.EMax(CL_idx(i), Drug_idx_B(i)) > T.Einf_B(i) )
            Mono.IC50(CL_idx(i), Drug_idx_B(i)) = T.IC50_B(i);
            Mono.EMax(CL_idx(i), Drug_idx_B(i)) = T.Einf_B(i);
            Mono.H(CL_idx(i), Drug_idx_B(i)) = T.H_B(i);        
            Mono.Max_C(CL_idx(i), Drug_idx_B(i)) = T.MAX_CONC_B(i);        
        end    
    end



    Mono.IC50(Mono.IC50==inf) = nan;
    Mono.EMax(Mono.EMax==inf) = nan;


    % TODO: Impute missing values based on Dual Layer method

    % Compute sensitivity score (AUC, atm)
    % TODO: CHECK ME!
    doses = [0,0.00001,0.00003,0.0001,0.0003,0.001
    0,0.00003,0.0001,0.0003,0.001,0.003
    0,0.0001,0.0003,0.001,0.003,0.01
    0,0.0003,0.001,0.003,0.01,0.03
    0,0.001,0.003,0.01,0.03,0.1
    0,0.003,0.01,0.03,0.1,0.3
    0,0.01,0.03,0.1,0.3,1
    0,0.03,0.1,0.3,1,3
    0,0.1,0.3,1,3,10
    0, 0.75, 2.5, 7.5, 25, 75];

    Mono.Drug_sensitivity = nan(size(annotations.cellLines, 1), size(annotations.drugs, 1));

    fun = @(a, IC50, EMax, h) 100 + ( (EMax - 100) ./ (1 + (IC50 ./ a).^h) );

    for i = 1:size(annotations.cellLines, 1)
        for j = 1:size(annotations.drugs, 1)
            if( ~isnan(Mono.IC50(i, j)) ) % && ~isnan(Mono.EMax(i, j)) && ~isnan(Mono.H(i, j)))
                current_dose = find(doses(:, end) == Mono.Max_C(i, j));
                if(isempty(current_dose))
                    fprintf('%d %d %e\n', i, j, Mono.Max_C(i, j));
                end
                Mono.Drug_sensitivity(i, j) = integral(@(x) fun(x, Mono.IC50(i, j), Mono.EMax(i, j), Mono.H(i, j)), doses(current_dose, 2), doses(current_dose, end));
            end
        end
    end

    Mono.IC50 = num2cell(Mono.IC50');
    Mono.H = num2cell(Mono.H');
    Mono.EMax = num2cell(Mono.EMax');
    Mono.Max_C = num2cell(Mono.Max_C');
    Mono.Drug_sensitivity = num2cell(Mono.Drug_sensitivity');
    
    %     [ii, jj] = find(~isnan(Mono.Drug_sensitivity));
    %     vv = Mono.Drug_sensitivity(sub2ind(size(Mono.Drug_sensitivity), ii, jj));
    %     zz = Modified_zscore(vv);

%         Mono.Drug_sensitivity = (Mono.Drug_sensitivity - nanmax(nonzeros(Mono.Drug_sensitivity))) ./ (nanmin(nonzeros(Mono.Drug_sensitivity)) - nanmax(nonzeros(Mono.Drug_sensitivity)));

end
