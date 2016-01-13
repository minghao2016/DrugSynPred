function [ C2C ] = C2C_expr( annotations, varargin )
    params = inputParser;
    params.addParamValue('expression_only'           , true, @(x) islogical(x) ); % should we use expression data only to compute C2C?
    params.addParamValue('selection_method'           , 'Elbow',@(x) ischar(x) ); % Cutting method for identifying selective genes using entropy
    params.addParamValue('SNF_K', 10      ,@(x) isscalar(x) & x > 0); %  Parameter K for the SNF (K-nearest neighbors only)
    params.addParamValue('SNF_T', 10      ,@(x) isscalar(x) & x > 0); %  Number of iterations for the SNF
    params.addParamValue('SNF_sigma', 0.5      ,@(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for computing Affinity matrix for SNF (width of RBF function)
    params.addParamValue('SNF_alpha', 1      ,@(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for computing Affinity matrix for SNF (width of RBF function)

    params.parse(varargin{:});
    par = params.Results;
    
    %% 1) C2C based on Dream microarray dataset
    [ Dream_expr ] = ReadExpr( 'Dream', annotations);

    [~,score] = pca(Dream_expr');
    C2C_Dream = corr(score');
    CellLine_distances{1} = 1 - C2C_Dream;
    
%     C2C_Dream(C2C_Dream < 0) = nan;
    
    
    %% 2) C2C based on CCLE microarray dataset
    [ CCLE_expr ] = ReadExpr( 'CCLE', annotations);

    [~,score] = pca(CCLE_expr');
    C2C_CCLE = corr(score');
    CellLine_distances{2} = 1 - C2C_CCLE;
    
%     C2C_CCLE(C2C_CCLE < 0) = nan;    

    
%     %% 3) C2C based on RNASeq
%     [ RNASeq_expr ] = ReadExpr( 'RNASeq_VSD', annotations);
% 
%     [~,score] = pca(RNASeq_expr');
%     C2C_RNASeq = corr(score');
%     CellLine_distances{3} = 1 - C2C_RNASeq;
%     
%     C2C_RNASeq(C2C_RNASeq < 0) = nan;    


        CellLine_distances{3} = ones(size(annotations.cellLines, 1));
        tissue_names = unique(annotations.cellLines.Disease_Area);
        for i = 1:numel(tissue_names)
            tissue_idx = find(strcmp( annotations.cellLines.Disease_Area, tissue_names{i}));
            CellLine_distances{3}(tissue_idx, tissue_idx) = 0;
        end
    %%
        W1 = affinityMatrix(CellLine_distances{1}, par.SNF_K, par.SNF_sigma);
        W2 = affinityMatrix(CellLine_distances{2}, par.SNF_K, par.SNF_sigma);

        par.SNF_K = 5;
        par.SNF_T = 10;
        par.alpha = 0.1;
        
        C2C = SNF({W1, W2}, par.SNF_K, par.SNF_T, par.SNF_alpha);    
%     %% 4) Combine
%     C2C = reshape(nanmean([C2C_Dream(:), C2C_CCLE(:), C2C_RNASeq(:)], 2), size(C2C_Dream, 1), size(C2C_Dream, 1));    
%     C2C(isnan(X)) = 0;
end

