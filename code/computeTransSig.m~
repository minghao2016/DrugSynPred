function     [ transcriptional_gene_signature, transc_class_gene_idx ] = computeTransSig( annotations, ACSN, varargin)
    params = inputParser;
    params.addParamValue('max_iter', 5, @(x) isscalar(x) & x > 0 & x <=1 ); % for double-propagation
    params.parse(varargin{:});
    par = params.Results;

    
    % TODO: Draft file! CHECK CHECK CHECK, to make sure we selected the best drugs to assay
    LINCS_ds = parse_gct('input/LINCS/final/LINCS_subset.gct');
    LINCS_genes = LINCS_ds.rdesc(:, 7);
    transc_class_gene_idx = cellfun(@(genes) find(ismember(LINCS_genes, genes)), ACSN.class_genes, 'UniformOutput', false) ;

    LINCS_celllines = LINCS_ds.cdesc(:, 1);
    LINCS_celllines(strcmp(LINCS_celllines, 'BT20')) = {'BT-20'};
    LINCS_celllines(strcmp(LINCS_celllines, 'HT29')) = {'HT-29'};    
    LINCS_celllines(strcmp(LINCS_celllines, 'MDAMB231')) = {'MDA-MB-231'};    
    LINCS_celllines(strcmp(LINCS_celllines, 'HS578T')) = {'Hs-578-T'};    
    
    LINCS_drugs = LINCS_ds.cdesc(:, 7);
    LINCS_expression_matrix = LINCS_ds.mat; % TODO: Should we column normalize to ensure constant transcriptional activity for each drug?   
    
    
    [~, cl_idx] = ismember(LINCS_celllines, annotations.cellLines.Sanger_Name);

    % TODO: Check to make sure all ID mappings are correct
    Dream2LINCS = readtable('./input/LINCS/final/preliminary_mapping.csv');
    
    transcriptional_gene_signature = cell(size(annotations.drugs, 1), size(annotations.cellLines, 1));
    for i = 1:size(LINCS_expression_matrix, 2)        
        drug_row_mask = ismember(Dream2LINCS.ID, LINCS_drugs{i});        
        transcriptional_gene_signature(drug_row_mask, cl_idx(i)) = {LINCS_expression_matrix(:, i)};
    end
    

%     %Propagate the expression per cell line and per drug to impute missing values
%     alpha = 0.1;
%     [D2D_prop_Exp,C2C_prop_Exp] = propagate_Expression(D2D, C2C , transcriptional_gene_signature, alpha);   
%     
%     
    % Cross-validate by masking 2 cell lines or 10 drugs
end

