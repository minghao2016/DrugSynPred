function [  ] = exportResults( annotations, Pairs, Pair_names, Confidence_mat,  outPath, varargin )
    fprintf('Exporting results to %s ... \n', outPath);
    
    params = inputParser;
    params.addParamValue('synergy_threshold', 30, @(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for random walk (larger alpha, deeper the length of random walks, i.e., alpha is the weight of topological similarity)
    params.parse(varargin{:});
    par = params.Results;

    sorted_CL = sort(annotations.cellLines.Sanger_Name);    
    fd_syn = fopen(fullfile(outPath, 'synergy_matrix.csv'), 'w');
    fd_conf = fopen(fullfile(outPath, 'confidence_matrix.csv'), 'w');
    for cIdx = 1:size(annotations.cellLines, 1)
        fprintf(fd_syn, ',%s', sorted_CL{cIdx});
        fprintf(fd_conf, ',%s', sorted_CL{cIdx});
    end
    fprintf(fd_syn, '\n');
    fprintf(fd_conf, '\n');    
    
    for pIdx = 1:size(Pairs, 1)
        fprintf(fd_syn, '%s.%s', Pair_names{pIdx, 1}, Pair_names{pIdx, 2});
        fprintf(fd_conf, '%s.%s', Pair_names{pIdx, 1}, Pair_names{pIdx, 2});
        for cIdx = 1:size(annotations.cellLines, 1)
            fprintf(fd_syn, ',%d', Confidence_mat(pIdx, cIdx) > par.synergy_threshold);
            fprintf(fd_conf, ',%f', Confidence_mat(pIdx, cIdx));
        end
        if(pIdx ~= size(Pairs, 1))
            fprintf(fd_syn, '\n');
            fprintf(fd_conf, '\n');    
        end
    end    
    fclose(fd_syn);
    fclose(fd_conf);

end

