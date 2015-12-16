function [ allMono ] = read_allMonoTherapy( annotations, inPath )

    if(~exist('input/preprocessed/allMono.mat', 'file'))
        fprintf('Importing allMono therapy data ...\n');
        % Read allMono therapy
        mono_files = dir(fullfile(inPath, 'ch*.csv'));        

        allMono.IC50 = inf(size(annotations.cellLines, 1), size(annotations.drugs, 1));
        allMono.EMax = inf(size(annotations.cellLines, 1), size(annotations.drugs, 1));
        allMono.H = nan(size(annotations.cellLines, 1), size(annotations.drugs, 1));
        allMono.Max_C = nan(size(annotations.cellLines, 1), size(annotations.drugs, 1));
        
        for f_id = 1:numel(mono_files)
            fname = mono_files(f_id).name;
            fprintf('\t%s\n', fname);
            T = readtable(fullfile(inPath, fname), 'TreatAsEmpty', 'NA', 'Format', '%s%s%s%f%f%f%f%f%f%f%f%f%f%s', 'Delimiter', ',');
            [~, CL_idx] = ismember(T.CELL_LINE, annotations.cellLines.Sanger_Name);
            [~, Drug_idx_A] = ismember(T.COMPOUND_A, annotations.drugs.ChallengeName);
            [~, Drug_idx_B] = ismember(T.COMPOUND_B, annotations.drugs.ChallengeName);
            for i = 1:size(T, 1)
                
                if( T.IC50_A(i) < allMono.IC50(CL_idx(i), Drug_idx_A(i)) && allMono.EMax(CL_idx(i), Drug_idx_A(i)) > T.Einf_A(i))
                    allMono.IC50(CL_idx(i), Drug_idx_A(i)) = T.IC50_A(i);
                    allMono.EMax(CL_idx(i), Drug_idx_A(i)) = T.Einf_A(i);
                    allMono.H(CL_idx(i), Drug_idx_A(i)) = T.H_A(i);        
                    allMono.Max_C(CL_idx(i), Drug_idx_A(i)) = T.MAX_CONC_A(i);        
                end
                if( T.IC50_B(i) < allMono.IC50(CL_idx(i), Drug_idx_B(i)) && allMono.EMax(CL_idx(i), Drug_idx_B(i)) > T.Einf_B(i) )
                    allMono.IC50(CL_idx(i), Drug_idx_B(i)) = T.IC50_B(i);
                    allMono.EMax(CL_idx(i), Drug_idx_B(i)) = T.Einf_B(i);
                    allMono.H(CL_idx(i), Drug_idx_B(i)) = T.H_B(i);        
                    allMono.Max_C(CL_idx(i), Drug_idx_B(i)) = T.MAX_CONC_B(i);        
                end    
            end
        end
            


        allMono.IC50(allMono.IC50==inf) = nan;
        allMono.EMax(allMono.EMax==inf) = nan;

        
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

        allMono.Drug_sensitivity = nan(size(annotations.cellLines, 1), size(annotations.drugs, 1));

        fun = @(a, IC50, EMax, h) 100 + ( (EMax - 100) ./ (1 + (IC50 ./ a).^h) );

        for i = 1:size(annotations.cellLines, 1)
            for j = 1:size(annotations.drugs, 1)
                if( ~isnan(allMono.IC50(i, j)) ) % && ~isnan(allMono.EMax(i, j)) && ~isnan(allMono.H(i, j)))
                    current_dose = find(doses(:, end) == allMono.Max_C(i, j));
                    if(isempty(current_dose))
                        fprintf('%d %d %e\n', i, j, allMono.Max_C(i, j));
                    end
                    allMono.Drug_sensitivity(i, j) = integral(@(x) fun(x, allMono.IC50(i, j), allMono.EMax(i, j), allMono.H(i, j)), doses(current_dose, 2), doses(current_dose, end));
                end
            end
        end


        %     [ii, jj] = find(~isnan(allMono.Drug_sensitivity));
        %     vv = allMono.Drug_sensitivity(sub2ind(size(allMono.Drug_sensitivity), ii, jj));
        %     zz = Modified_zscore(vv);

%         allMono.Drug_sensitivity = (allMono.Drug_sensitivity - nanmax(nonzeros(allMono.Drug_sensitivity))) ./ (nanmin(nonzeros(allMono.Drug_sensitivity)) - nanmax(nonzeros(allMono.Drug_sensitivity)));

        allMono.IC50 = num2cell(allMono.IC50');
        allMono.H = num2cell(allMono.H');
        allMono.EMax = num2cell(allMono.EMax');
        allMono.Max_C = num2cell(allMono.Max_C');
        allMono.Drug_sensitivity = num2cell(allMono.Drug_sensitivity');

        save('input/preprocessed/allMono.mat', 'allMono');
    else
        fprintf('Loading allMono therapy data ...\n');
        load('input/preprocessed/allMono.mat');
    end
end

