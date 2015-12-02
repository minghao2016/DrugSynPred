function [ C2C, Expressed_genes ] = Construct_C2C( annotations, varargin )
    params = inputParser;
    params.addParamValue('expression_only'           , true, @(x) islogical(x) ); % should we use expression data only to compute C2C?
    params.addParamValue('selection_method'           , 'Elbow',@(x) ischar(x) ); % Cutting method for identifying selective genes using entropy
    params.addParamValue('SNF_K', 10      ,@(x) isscalar(x) & x > 0); %  Parameter K for the SNF (K-nearest neighbors only)
    params.addParamValue('SNF_T', 10      ,@(x) isscalar(x) & x > 0); %  Number of iterations for the SNF
    params.addParamValue('SNF_alpha', 0.5      ,@(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for computing Affinity matrix for SNF (width of RBF function)
        
    params.parse(varargin{:});
    par = params.Results;

    
    m = size(annotations.cellLines, 1); % # of cell lines
    n = size(annotations.drugs, 1); % # of drugs
    
    CellLine_distances = cell(3, 1);

if(~exist('input/preprocessed/C2C.mat', 'file'))    
    %% Tissue-based similarity
    CellLine_distances{1} = ones(size(annotations.cellLines, 1));
    tissue_names = unique(annotations.cellLines.Tissue__General_);
    for i = 1:numel(tissue_names)
        tissue_idx = find(strcmp( annotations.cellLines.Tissue__General_, tissue_names{i}));
        CellLine_distances{1}(tissue_idx, tissue_idx) = 0;
    end
    CellLine_distances{1} = CellLine_distances{1};
    
    %% Process Expression datasets
    % Co-expr matrix from CCLE
    fd = fopen('input/CCLE/CCLE_Expression_Entrez_2012-09-29.gct', 'r');
    if(fd == -1)
        error('Cannot open file input/CCLE/CCLE_Expression_Entrez_2012-09-29.gct');
    else        
        header = fscanf(fd, '%[^\r\n]\r\n', 1);
        cols = textscan(header, '%s', 'Whitespace', '', 'Delimiter', '\t');
        colNames = cols{1};

        FormatString = sprintf('%%s%%s%s%%[^\\r\\n]', repmat('%f',1, numel(colNames)-3));        
        table = textscan(fd, FormatString, 'Whitespace', '', 'Delimiter', '\t');
    end
    fclose(fd);
    gene_names = table{2};    
    
    CCLE_mask = ismember(colNames, annotations.cellLines.CCLE_Name);
    table = table(CCLE_mask);
    colNames = colNames(CCLE_mask);
    
    CCLE_expr = zeros(numel(table{1}), size(annotations.cellLines, 1));    
    [~, selected_cols] = ismember(colNames, annotations.cellLines.CCLE_Name);
    CCLE_expr(:, selected_cols) = cell2mat(table);
    CCLE_expressed = arrayfun(@(col) gene_names((CCLE_expr(:, col) > 6)), 1:size(annotations.cellLines, 1), 'UniformOutput', false);
    
    P = normalize(CCLE_expr, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
    I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
    H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
    sorted_H = sort(H);
    sigSelective_no = cut(H, 'selection_method', par.selection_method, 'selection_half', 'bottom');
    entropy_cutoff = sorted_H(sigSelective_no);
    
    expression_mat = CCLE_expr(H < entropy_cutoff, :);
    
    [CC, CC_pval] = corr(expression_mat);   
    CC_pval(CC < 0) = 1;
    CC_pval(1e-10 < CC_pval) = 1;
    CC_pval(CC_pval == 0) = min(nonzeros(CC_pval));    
    S = -log10(CC_pval);

    [ii, jj, vv] = find(S);
    D = full(sparse(ii, jj, (max(vv) - vv) ./ (max(vv) - min(vv)), size(S, 1), size(S, 1)));
    D = D - diag(diag(D));  
    D(isnan(D)) = 1;
    co_exp{1} = D; 

    

    % Co-expr matrix from Dream
    [expr_table, cellLine_names, rowNames] = my_tblread('input/Dream/molecular/gex.csv', ',');
    cellLine_names = cellfun(@(x) x(2:end-1), cellLine_names, 'UniformOutput', false);

    Dream_expr = zeros(size(expr_table, 1), size(annotations.cellLines, 1));    
    [~, selected_cols] = ismember(cellLine_names, annotations.cellLines.Sanger_Name);
    Dream_expr(:, selected_cols) = expr_table;
    Dream_expressed = arrayfun(@(col) rowNames((Dream_expr(:, col) > 5)), 1:size(annotations.cellLines, 1), 'UniformOutput', false);
    Expressed_genes = cellfun(@(x, y) union(x, y), Dream_expressed, CCLE_expressed, 'UniformOutput', false);

    P = normalize(Dream_expr, 'dim', 2, 'pnorm', 1); % Relative expression of each gene in different tissues
    I = spfun(@(x) -log2(x), P); % Information content of the relative expression of each gene   
    H = sum(P .* I, 2); % Shannon entropy of the relative expression -- overall tissue-specificity of genes. It has unit of "bits" and is between 0 -> completely selective) and log2(k) -> Uniform/HK gene    
    sorted_H = sort(H);
    sigSelective_no = cut(H, 'selection_method', par.selection_method, 'selection_half', 'bottom');
    entropy_cutoff = sorted_H(sigSelective_no);
    
    expression_mat = Dream_expr(H < entropy_cutoff, :);
    
    [CC, CC_pval] = corr(expression_mat);   
    CC_pval(CC < 0) = 1;
    CC_pval(1e-10 < CC_pval) = 1;
    CC_pval(CC_pval == 0) = min(nonzeros(CC_pval));    
    S = -log10(CC_pval);

    [ii, jj, vv] = find(S);
    D = full(sparse(ii, jj, (max(vv) - vv) ./ (max(vv) - min(vv)), size(S, 1), size(S, 1)));
    D = D - diag(diag(D));  
    D(isnan(D)) = 1;
    co_exp{2} = D; 

        
    
    M = min(co_exp{1}, co_exp{2});
    sigma = std(nonzeros(tril(M)));
    [ii, jj, vv] = find(tril(M));
    vv = arrayfun(@(x) exp(x / (2*sigma)), vv);
    vv = vv / max(vv);
    combined_D = sparse(ii, jj, vv, m, m);
    combined_D = max(combined_D, combined_D');
    CellLine_distances{2} = combined_D;
    
    if(par.expression_only == true)
        C2C = 1 - CellLine_distances{2};

%         breast_subtype = {'BT-20','Basal';'BT-474','Luminal';'BT-549','Claudin-low';'HCC1143','Basal';'HCC1187','Basal';'HCC1395','Claudin-low';'HCC1419','Luminal';'HCC1428','Luminal';'HCC1806','Basal';'HCC1937','Basal';'HCC1954','Basal';'HCC38','Claudin-low';'HCC70','Basal';'Hs-578-T','Claudin-low';'MCF7','Luminal';'MDA-MB-157','Claudin-low';'MDA-MB-231','Claudin-low';'MDA-MB-361','Luminal';'MDA-MB-415','Luminal';'MDA-MB-436','Claudin-low';'MDA-MB-453','Luminal';'MDA-MB-468','Basal';'T47D','Luminal'}; % From Sichani or Heiser?
%         Breast = C2C(1:35, 1:35);            
%         avg = mean(nonzeros(tril(Breast)))
%         types = unique(breast_subtype(:, 2));
%         for i = 1:numel(types)
%             class_CLs = breast_subtype(strcmp(types{i}, breast_subtype(:, 2)), 1);            
%             [~, idx] = ismember(class_CLs, annotations.cellLines.Sanger_Name);
%             class_avg = mean(nonzeros(tril(Breast(idx, idx))))
%         end
        
        save('input/preprocessed/C2C.mat', 'expr_table', 'H', 'C2C', 'Expressed_genes');        
        return;
    end
    

%     clustergram(CellLine_distances{2}, 'RowLabels', annotations.cellLines.CCLE_Name, 'ColumnLabels', annotations.cellLines.CCLE_Name, 'Linkage', 'average', 'ColorMap', colormap(flipud(redgreencmap())), 'OPTIMALLEAFORDER', true)    

    %% Process mutations   
    tic; [MAF] = my_dlmread('input/Dream/molecular/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.txt', '\t'); toc

    Census_genes = readtable('input/Dream/molecular/Cancer_Census_genes.txt', 'Delimiter', '\t');    
    Census_mask = ismember(MAF(:, 1), Census_genes.HGNCSymbol);    
    MAF(~Census_mask, :) = [];

    Census_map = zeros(numel(Census_genes), m);
    [~, gene_idx] = ismember(MAF(:, 1), Census_genes.HGNCSymbol);    
    [~, cl_idx] = ismember(MAF(:, 16), annotations.cellLines.CCLE_Name);
    
    matched_cols = setdiff(unique(cl_idx), 0);
    for i = 1:numel(matched_cols)
        col = matched_cols(i);
        mutated_genes = gene_idx(cl_idx == col, :);
        signature = accumarray(mutated_genes, 1);
        [ii, ~, vv] = find(signature);
        Census_map(ii, col) = vv;        
    end
 
    D = squareform(pdist(double(logical(Census_map))', 'Hamming'));
    [ii, jj, vv] = find(D);
    D = sparse(ii, jj, (vv - min(vv)) ./ (max(vv) - min(vv)), size(D, 1), size(D, 1));
    D(:, setdiff(1:m, matched_cols)) = 1;
    D(setdiff(1:m, matched_cols), :) = 1;   
    CellLine_distances{3} = D;
%     clustergram(CellLine_distances{3}, 'RowLabels', annotations.cellLines.CCLE_Name, 'ColumnLabels', annotations.cellLines.CCLE_Name, 'Linkage', 'average', 'OPTIMALLEAFORDER', true)    

    %% Combine networks
    % Run SNF

    % TODO: 1) Should we normalize distances first?, 2) Should we
    % incorporate mutations, CNV, mutation...? Which one of them?, 3) What
    % are the optimal parameters for SNF? 
%     W1 = CellLine_distances{1};
%     W2 = CellLine_distances{2};
%     W3 = CellLine_distances{3};
    
    W1 = affinityMatrix(CellLine_distances{1}, par.SNF_K, par.SNF_alpha);
    W2 = affinityMatrix(CellLine_distances{2}, par.SNF_K, par.SNF_alpha);
    W3 = affinityMatrix(CellLine_distances{2}, par.SNF_K, par.SNF_alpha);
    
    C2C = SNF({W2, W3}, par.SNF_K, par.SNF_alpha);

    C2C = (C2C - min(nonzeros(C2C))) ./ (max(nonzeros(C2C)) - min(nonzeros(C2C)));
    C2C = C2C - diag(diag(C2C));
%     C2C = affinityMatrix(C2C, par.SNF_K, par.SNF_alpha); % TODO: Should we use affinityMatrix at the end?
%     clustergram(C2C, 'RowLabels', annotations.cellLines.CCLE_Name, 'ColumnLabels', annotations.cellLines.CCLE_Name, 'Linkage', 'average', 'OPTIMALLEAFORDER', true)    


    save('input/preprocessed/C2C.mat', 'expr_table', 'H', 'MAF', 'Census_genes', 'CellLine_distances', 'C2C', 'Expressed_genes');

else
    load('input/preprocessed/C2C.mat');
end
end

