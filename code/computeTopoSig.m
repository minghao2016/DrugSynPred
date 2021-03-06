function [ topological_gene_signature, topo_class_gene_idx ] = computeTopoSig( annotations, ACSN, NodeWeights, varargin)
    params = inputParser;
    params.addParamValue('alpha', 0.85, @(x) isscalar(x) & x > 0 & x <=1 ); % RWR: for constructing TS network
    params.addParamValue('beta', 0.85, @(x) isscalar(x) & x > 0 & x <=1 ); % RWR: for computing topological signatures starting from drug targets
        
    params.parse(varargin{:});
    par = params.Results;
    
    fname = sprintf('input/preprocessed/topoSigs_beta=%.2f.mat', par.beta);
    if(~exist(fname, 'file'))
        fprintf('Computing topological gene signatures ...\n');
        n = size(ACSN.A, 1);
        W = ACSN.A;       

        topological_gene_signature = cell(size(annotations.drugs, 1), size(annotations.cellLines, 1));

        for j = 1:size(annotations.cellLines)
            fprintf('\t%d- Cell Line %s ...\n', j, annotations.cellLines.Sanger_Name{j});

            % Construct cell type-specific network
            fprintf('\tComputing cell-type specific network ...\n');
            v_score = NodeWeights(:, j);       
            smoothed_penalty = (v_score * size(annotations.cellLines, 1)).^2;
            CL_A = bsxfun(@times, W, smoothed_penalty'); % Penalize rows (sources)            
            CL_A = bsxfun(@times, CL_A, smoothed_penalty); % Penalize columns (destinations)            

            
            % Perform random walk on cell-type-scpecific network starting from
            % the targets of each drug
            fprintf('\tPrecomputing random-walk matrix ...\n');
            P = CL_A'*spdiags(spfun(@(x) 1./x, sum(CL_A, 2)), 0, n, n);
            e_T = ones(1, n);    
            d_T = e_T - e_T*P;
            P = P + diag(d_T);
            Q = (1-par.beta).*inv(eye(n) - par.beta*P);

            for i = 1:size(annotations.drugs, 1)
                fprintf('\t\t%d- Drug %s ...\n', i, annotations.drugs.ChallengeName{i});
                primary_targets = annotations.drugs.Target{i};
                src_nodes = find(cellfun(@(x) nnz(ismember(primary_targets, x)), ACSN.vertex_genes));
                if(numel(src_nodes) == 0)
                    continue;
                end
                e_src = sparse(src_nodes, 1, 1, n, 1); e_src = e_src ./ sum(e_src);
                topological_gene_signature{i, j} = Q*e_src;
            end
        end
        topo_class_gene_idx = cellfun(@(genes) find(cellfun(@(v) nnz(ismember(v, genes)), ACSN.vertex_genes)), ACSN.class_genes, 'UniformOutput', false) ; % For topological
        save(fname, 'topological_gene_signature', 'topo_class_gene_idx');
    else
        fprintf('Loading topological gene signatures ...\n');
        load(fname, 'topological_gene_signature', 'topo_class_gene_idx');
    end
end

