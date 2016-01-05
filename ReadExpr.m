function [ Expr_mat, gene_names, expr_threshold, target_idx ] = ReadExpr( name, annotations)

    switch(name)
        case 'RNASeq_VSD'
            [T] = my_tdfread('input/CL_RNASeq/E-MTAB-2706.sdrf.txt'); 

            fd = fopen('input/CL_RNASeq/Dream2RNASeq.csv', 'r');
            C = textscan(fd, '%s %s', 'Delimiter', ',');
            fclose(fd);

            tic;
            [VSD, samples, genes] = my_tblread('input/CL_RNASeq/140625_Klijn_RPKM_coding.txt');
            toc

            [~, sample_idx] = ismember(samples, T{1});
            sample_names = T{7}(sample_idx);
            selection_mask = ismember(sample_names, C{2});
            sample_names(~selection_mask) = [];
            VSD(:, ~selection_mask) = [];

            Expr_mat = zeros(size(VSD, 1), size(annotations.cellLines, 1));
            [~, col_idx] = ismember(sample_names, C{2});
            Expr_mat(:, col_idx) = VSD;            
            
            [G] = my_tdfread('input/CL_RNASeq/Entrez2HUGO.txt');
            [row_mask] = ismember(genes, G{2});
            Expr_mat(~row_mask, :) = [];
            genes(~row_mask) = [];
            [~, idx] = ismember(genes, G{2});
            gene_names = G{1}(idx);
                
            target_idx = arrayfun(@(x) find(ismember(gene_names, annotations.drugs.Target{x})), 1:size(annotations.drugs, 1), 'UniformOutput', false);
            
            expr_threshold = 500;
            
        case 'RNASeq_RPKM'
            [T] = my_tdfread('input/CL_RNASeq/E-MTAB-2706.sdrf.txt'); 

            fd = fopen('input/CL_RNASeq/Dream2RNASeq.csv', 'r');
            C = textscan(fd, '%s %s', 'Delimiter', ',');
            fclose(fd);

            tic;
            [RPKM, samples, genes] = my_tblread('input/CL_RNASeq/140625_Klijn_VSD_coding.txt');
            toc

            [~, sample_idx] = ismember(samples, T{1});
            sample_names = T{7}(sample_idx);
            selection_mask = ismember(sample_names, C{2});
            sample_names(~selection_mask) = [];
            RPKM(:, ~selection_mask) = [];

            Expr_mat = zeros(size(RPKM, 1), size(annotations.cellLines, 1));
            [~, col_idx] = ismember(sample_names, C{2});
            Expr_mat(:, col_idx) = RPKM;      

            [G] = my_tdfread('input/CL_RNASeq/Entrez2HUGO.txt');
            [row_mask] = ismember(genes, G{2});
            Expr_mat(~row_mask, :) = [];
            genes(~row_mask) = [];
            [~, idx] = ismember(genes, G{2});
            gene_names = G{1}(idx);
                
            target_idx = arrayfun(@(x) find(ismember(gene_names, annotations.drugs.Target{x})), 1:size(annotations.drugs, 1), 'UniformOutput', false);

            expr_threshold = 1;
            
        case 'Dream'
            [expr_table, cellLine_names, gene_names] = my_tblread('input/Dream/molecular/gex.csv', ',');
            cellLine_names = cellfun(@(x) x(2:end-1), cellLine_names, 'UniformOutput', false);
            gene_names = cellfun(@(x) x(2:end-1), gene_names, 'UniformOutput', false);

            Expr_mat = zeros(size(expr_table, 1), size(annotations.cellLines, 1));    
            [~, selected_cols] = ismember(cellLine_names, annotations.cellLines.Sanger_Name);

            Expr_mat(:, selected_cols) = expr_table; %zscore(expr_table, 0, 2);

            target_idx = arrayfun(@(x) find(ismember(gene_names, annotations.drugs.Target{x})), 1:size(annotations.drugs, 1), 'UniformOutput', false);
            
            expr_threshold = 5.5;

        case 'CCLE'
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

            Expr_mat = zeros(numel(table{1}), size(annotations.cellLines, 1));    
            [~, selected_cols] = ismember(colNames, annotations.cellLines.CCLE_Name);
            Expr_mat(:, selected_cols) = cell2mat(table);  

            expr_threshold = 12;            
    end
            
end

