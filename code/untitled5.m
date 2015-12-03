function     [ transcriptional_gene_signature ] = computeTanscSig( annotations, ACSN, Expressed_genes, varargin)
    params = inputParser;
    params.addParamValue('synergy_threshold', 30, @(x) isscalar(x) & x > 0 & x <=1 ); % Alpha parameter for random walk (larger alpha, deeper the length of random walks, i.e., alpha is the weight of topological similarity)
    params.parse(varargin{:});
    par = params.Results;

end

