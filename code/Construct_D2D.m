function [ D2D, coreD2D, targetD2D, topological_signatures ] = Construct_D2D(annotations, interactome, varargin)
    params = inputParser;
    params.addParamValue('alpha', 0.9, @(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for random walk (larger alpha, deeper the length of random walks, i.e., alpha is the weight of topological similarity)
    params.addParamValue('lambda', 0.5, @(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for random walk (larger alpha, deeper the length of random walks, i.e., alpha is the weight of topological similarity)
        
    params.parse(varargin{:});
    par = params.Results;
    
    if(~exist('input/preprocessed/D2D.mat', 'file'))
        fprintf('Constructing D2D ...\n');
    % Exand targets using WinDTome
%         m = size(annotations.drugs, 1);
    %     % Read drug targets
    %     if(~exist('input/Drug2Target/WinDTome/WinDTome.mat', 'file'))
    %         tic; 
    %         WinDTome = readtable('input/Drug2Target/WinDTome/WinDTome.txt', 'Delimiter', '\t', 'Format', '%s %s %d %s %s %s %s %s %s %s %d %d'); 
    % 
    %     %     % based on Drug name
    %     %     [WinDTome_drugs, ~, ia] = unique(WinDTome.Drug_ID);
    %     %     WinDTome_targets = arrayfun(@(drug_id) unique(WinDTome.Target_Gene_Symbol(ia == drug_id)), 1:numel(WinDTome_drugs), 'UniformOutput', false)';
    % 
    %         % Based on DrugBank ID
    %         WinDTome(strcmp(WinDTome.Source_DrugBank, ''), :) = [];
    % 
    %         [WinDTome_drugs, ~, ia] = unique(WinDTome.Source_DrugBank);
    %         WinDTome_targets = arrayfun(@(drug_id) unique(WinDTome.Target_Gene_Symbol(ia == drug_id)), 1:numel(WinDTome_drugs), 'UniformOutput', false)';
    % 
    % 
    %     %     drug_homologs = cell(size(annotations.drugs, 1), 1);
    %     %     for i = 1:size(annotations.drugs, 1)
    %     %         targets = annotations.drugs.Target{i};
    %     % 
    %     %     %     target_overlap = cellfun(@(D) numel(intersect(D, targets)) / numel(union(D, targets)), WinDTome_targets);
    %     %     %     target_overlap = cellfun(@(D) numel(intersect(D, targets)), WinDTome_targets);
    %     %         tic; target_overlap = arrayfun(@(d) nnz(ismember(annotations.drugs.Target{i}, WinDTome_targets{d})), 1:numel(WinDTome_drugs)); toc
    %     %         [~, drug_row, drug_overlap] = find(target_overlap);
    %     %         if(numel(drug_row) == 0)
    %     %             drug_homologs{i} = {};
    %     %         else
    %     %             homologous_drug_IDs = WinDTome_drugs(drug_row);
    %     %             [~, perm] = sort(zscore(drug_overlap), 'descend');
    %     %             homologous_drug_IDs = homologous_drug_IDs(perm);
    %     %             drug_homologs{i} = homologous_drug_IDs;
    %     %         end
    %     %     end
    %     %     toc
    %         save('input/Drug2Target/WinDTome/WinDTome.mat', 'WinDTome', 'WinDTome_drugs', 'WinDTome_targets');
    %     else
    %         load('input/Drug2Target/WinDTome/WinDTome.mat'); 
    %     end
    

    % Read core D2D from STITCH network
        stitch = load('input/STITCH/drug_comb_mapped');
        coreD2D = sparse(stitch(:,1), stitch(:,2), stitch(:, 3), 119, 119);
        coreD2D = max(coreD2D, coreD2D');
        coreD2D = full(coreD2D ./ max(coreD2D(:)));
        coreD2D(~coreD2D) = nan;

    % Compute RWR over interactome from targets to identify topological signature of each drug
        n = size(interactome.A, 1);
        P = interactome.A'*spdiags(spfun(@(x) 1./x, sum(interactome.A, 2)), 0, n, n);

        % Handling the dangling nodes ...
        e_T = ones(1, n);    
        d_T = e_T - e_T*P;
        P = P + diag(d_T);
        Q = (1-par.alpha).*inv(eye(n) - par.alpha*P);

        topological_signatures = zeros(n, size(annotations.drugs, 1));
        for i = 1:size(annotations.drugs, 1)
            primary_targets = annotations.drugs.Target{i};
            [~, src_nodes] = ismember(primary_targets, interactome.vertex_genes);
            src_nodes(~src_nodes) = [];
            if(numel(src_nodes) == 0)
                continue;
            end
            e_src = sparse(src_nodes, 1, 1, n, 1); e_src = e_src ./ sum(e_src);
            topological_signatures(:, i) = Q*e_src;
        end
        targetD2D = partialcorr(topological_signatures, mean(topological_signatures, 2));
        pos_targetD2D = targetD2D;
        pos_targetD2D(pos_targetD2D < 0) = nan;
        pos_targetD2D = pos_targetD2D - diag(diag(pos_targetD2D));
        
        
%         targetD2D = topological_signatures'*topological_signatures;            
%         targetD2D = targetD2D ./ max(targetD2D(:));
%         X=log10(topological_signatures);
%         X(isinf(X)) = inf;
%         X = X - min(nonzeros(X));
%         X(isinf(X)) = 0;
%         X(sum(X, 2) == 0, :) = [];
%         Y = 1+corr(X);
%         [ii, jj, vv] = find(tril(targetD2D));
%         sigma = std(vv);
%         vv = arrayfun(@(x) exp(x / (6*sigma)), vv);
%         vv = vv / max(vv);
%         targetD2D = sparse(ii, jj, vv, m, m);
%         targetD2D = max(targetD2D, targetD2D');
    %     
    %     X = Modified_zscore();

%         D2D = par.lambda*coreD2D + (1-par.lambda)*targetD2D;  
        D2D = max(coreD2D, pos_targetD2D);

        save('input/preprocessed/D2D.mat', 'D2D', 'coreD2D', 'targetD2D', 'topological_signatures');
    else
        fprintf('Loading D2D ...\n');
        load('input/preprocessed/D2D.mat');
    end
end

