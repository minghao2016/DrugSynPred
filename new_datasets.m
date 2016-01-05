    clear
    addpath(genpath('code'));
    warning('off','all');    

    annotations.cellLines = readtable('input/Dream/molecular/cell_info.csv', 'Delimiter', ',');
    annotations.drugs = readtable('input/Dream/synergy/Drugs_final.txt', 'Delimiter', '\t');
    annotations.drugs.Target = cellfun(@(targets) strsplit(targets, ';'), annotations.drugs.Target, 'UniformOutput', false);

    [~, CL_perm] = sort(annotations.cellLines.Tissue__General_);
    annotations.cellLines = annotations.cellLines(CL_perm, :);
    
%%    
    CTRC_ds = parse_gct('input/Cancer Therapeutics Response Portal/Cancer_Therapeutics_Response_Cluster.gct');


%%
    [T, Header] = my_tdfread('input/CL_RNASeq/E-MTAB-2706.sdrf.txt'); 
    
    fd = fopen('input/CL_RNASeq/Dream2RNASeq.csv', 'r');
    C = textscan(fd, '%s %s', 'Delimiter', ',');
    fclose(fd);

    tic;
    [VSD, samples, genes] = my_tblread('input/CL_RNASeq/140625_Klijn_RPKM_coding.txt');
%     [RPKM, samples, genes] = my_tblread('input/CL_RNASeq/140625_Klijn_VSD_coding.txt');
    toc
    
    [~, sample_idx] = ismember(samples, T{1});
    sample_names = T{7}(sample_idx);
    selection_mask = ismember(sample_names, C{2});
    sample_names(~selection_mask) = [];
    VSD(:, ~selection_mask) = [];
    
    RNASeq_VSD_expr = zeros(size(VSD, 1), size(annotations.cellLines, 1));
    [~, col_idx] = ismember(sample_names, C{2});
    RNASeq_VSD_expr(:, col_idx) = VSD;

    [G] = my_tdfread('input/CL_RNASeq/Entrez2HUGO.txt');
    [row_mask, row_idx] = ismember(genes, G{2});
    
%%
        [expr_table, cellLine_names, Dream_gene_names] = my_tblread('input/Dream/molecular/gex.csv', ',');
        cellLine_names = cellfun(@(x) x(2:end-1), cellLine_names, 'UniformOutput', false);
        Dream_gene_names = cellfun(@(x) x(2:end-1), Dream_gene_names, 'UniformOutput', false);

        Dream_expr = zeros(size(expr_table, 1), size(annotations.cellLines, 1));    
        [~, selected_cols] = ismember(cellLine_names, annotations.cellLines.Sanger_Name);

        Dream_expr(:, selected_cols) = expr_table; %zscore(expr_table, 0, 2);

        drug_target_expr_indices = arrayfun(@(x) find(ismember(Dream_gene_names, annotations.drugs.Target{x})), 1:size(annotations.drugs, 1), 'UniformOutput', false);

    
    
    
    
    