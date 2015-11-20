function newkernels = combineKernels(geneList, kernels)
%function varargout = combineKernels(geneList, kernels)
% geneList is a unique list of genes which is a superset of all genes that
% are in kernels (i.e. all genes in kernels must be present in geneList)
% kernels is a cell structure of kernels whereby each entry of kernels has
% the fields kernels{1}.data, kernels{1}.rowlabels and
% kernels{1}.collabels; 
% newkernels is a cell structure of kernels with the same number of kernels
% as in kernels, now each kernel has the same number of genes ordered as in
% geneList

NN = length(geneList);
for ii = 1:length(kernels)
    K = kernels{ii};
    ix = lookup(K.rowlabels, geneList);
    if any(isnan(ix))
        error('Unmatched gene name');
    end
    Knew.rowlabels = geneList;
    Knew.collabels = geneList;
    [rr, cc, ss] = find(K.data);
    Knew.data = sparse(ix(rr), ix(cc), ss, NN, NN);
    newkernels{ii} = Knew;
end
