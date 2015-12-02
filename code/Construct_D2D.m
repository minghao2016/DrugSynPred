function [ D2D ] = Construct_D2D(ACSN, annotations, varargin)
    params = inputParser;
    params.addParamValue('alpha', 0.85, @(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for random walk (larger alpha, deeper the length of random walks, i.e., alpha is the weight of topological similarity)
        
    params.parse(varargin{:});
    par = params.Results;
    
% Read drug targets
if(~exist('input/Drug2Target/WinDTome/WinDTome.mat', 'file'))
    tic; 
    WinDTome = readtable('input/Drug2Target/WinDTome/WinDTome.txt', 'Delimiter', '\t', 'Format', '%s %s %d %s %s %s %s %s %s %s %d %d'); 
    
%     % based on Drug name
%     [WinDTome_drugs, ~, ia] = unique(WinDTome.Drug_ID);
%     WinDTome_targets = arrayfun(@(drug_id) unique(WinDTome.Target_Gene_Symbol(ia == drug_id)), 1:numel(WinDTome_drugs), 'UniformOutput', false)';

    % Based on DrugBank ID
    WinDTome(strcmp(WinDTome.Source_DrugBank, ''), :) = [];
    
    [WinDTome_drugs, ~, ia] = unique(WinDTome.Source_DrugBank);
    WinDTome_targets = arrayfun(@(drug_id) unique(WinDTome.Target_Gene_Symbol(ia == drug_id)), 1:numel(WinDTome_drugs), 'UniformOutput', false)';
    

%     drug_homologs = cell(size(annotations.drugs, 1), 1);
%     for i = 1:size(annotations.drugs, 1)
%         targets = annotations.drugs.Target{i};
% 
%     %     target_overlap = cellfun(@(D) numel(intersect(D, targets)) / numel(union(D, targets)), WinDTome_targets);
%     %     target_overlap = cellfun(@(D) numel(intersect(D, targets)), WinDTome_targets);
%         tic; target_overlap = arrayfun(@(d) nnz(ismember(annotations.drugs.Target{i}, WinDTome_targets{d})), 1:numel(WinDTome_drugs)); toc
%         [~, drug_row, drug_overlap] = find(target_overlap);
%         if(numel(drug_row) == 0)
%             drug_homologs{i} = {};
%         else
%             homologous_drug_IDs = WinDTome_drugs(drug_row);
%             [~, perm] = sort(zscore(drug_overlap), 'descend');
%             homologous_drug_IDs = homologous_drug_IDs(perm);
%             drug_homologs{i} = homologous_drug_IDs;
%         end
%     end
%     toc
    save('input/Drug2Target/WinDTome/WinDTome.mat', 'WinDTome', 'WinDTome_drugs', 'WinDTome_targets');
else
    load('input/Drug2Target/WinDTome/WinDTome.mat'); 
end

    
% Exand targets using WinDTome



% Compute RWR over ACSN from targets to identify topological signature of each drug
    fprintf('Computing random walk matrix ...\n');
    tic
    n = size(ACSN.A, 1);
    P = ACSN.A'*spdiags(spfun(@(x) 1./x, sum(ACSN.A, 2)), 0, n, n);
    
    % Handling the dangling nodes ...
    e_T = ones(1, n);    
    d_T = e_T - e_T*P;
    P = P + diag(d_T);
    Q = (1-par.alpha).*inv(eye(n) - par.alpha*P);
    toc
    
    topological_signatures = zeros(n, size(annotations.drugs, 1));
    for i = 1:size(annotations.drugs, 1)
        primary_targets = annotations.drugs.Target{i};
        src_nodes = find(cellfun(@(x) nnz(ismember(primary_targets, x)), ACSN.vertex_genes));
        if(numel(src_nodes) == 0)
            continue;
        end
        e_src = sparse(src_nodes, 1, 1, n, 1); e_src = e_src ./ sum(e_src);
        topological_signatures(:, i) = Q*e_src;
    end
    topological_signatures = topological_signatures / norm(topological_signatures);
    D2D = topological_signatures'*topological_signatures;
end

