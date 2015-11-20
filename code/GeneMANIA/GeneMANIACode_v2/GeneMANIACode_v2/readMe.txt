This folder contains the basic GeneMANIA code for predicting a function class using networks of gene-gene association data. 

*** To use the GeneMANIA code the following inputs are needed: 

1. A cell structure where each entry is an association networks:
   i.e. each entry in the cell structure is an NxN matrix W (where N is the number of genes), and 
   entry W(i,j) is the interactiton score between gene i and gene j. W must be symmetric all the entries in 
   W must be greater than or equal to 0. 

2. label vector of length Nx1 where N is the number of genes in the association networks. Entries of N are {-1, 0, +1}
   where genes that are positive for a given function class are +1, genes that are negative examples are 0, and unannotated
   genes are 0. 

*********
runGeneMANIA.m is a sample script which shows how to create association 
networks from text files and run GeneMANIA on a set of networks. 
**********

List of MATLAB m files for use with GeneMANIA: 

1. calcROCarea.m
    is required to calculate the ROC values for cross-validation
2. conjGrad.m
    this is the conjugate gradient code for solving a linear system of equations. 
3. findKernelWeights.m
    Given a MATLAB structure of association networks and a label network, this function finds a weight for each network and
    returns [alpha,indices] where alpha is the network weights, indices are the index of the association networks in the 
    input structure that are assigned weights in alpha (bias can be ignored). 
4. getScoreVectorCG.m
    This function computes the Laplacian of the composite association network and calls conjugate gradient. 
5. listResolve.m and lookup.m 
    These are auxilary files that are used by combineKernels.m
6. predictClassesCG.m
    Given a MATLAB structure with association networks, a label vector, returs a real-valued score
    vector that can be thresholded appropriately (no thresholding is applied by this function), K: a list of indices of 
    the kernels which received nonzero weight, and WTS, the corresponding weights assigned to the kernels 
    indexed in K (in the same order as K).
7. predictClassesCV.m
    Returns AREAS, a vector of length NFOLDS containing the ROC area on each step
    of cross-validation. Optionally returns R, a cell array of length NFOLDS, 
    containing the score vectors for each cross-validation step, and B, which
    returns a score vector of leave-out scores for each gene (i.e., the score
    that gene got when it was one of the left out genes).


List of auxilary files for preparing association networks:
1. makeAssociationKernel.m
    Given an MxN data matrix with N genes and M features and an integer K, 
    it returns an NxN matrix G where G(I,J) = 0 if J is not within the K closest neighbours
    otherwise, G(I,J) = the Pearson correlation of X(:,I) and X(:,J)
2. makeBinaryKernel.m
    same as above but more appropriate to use with binary data. 
    Here, G(I,J) = "Fisher kernel" of X(:,I) and X(:,J)
    in this case, the Fisher kernel is the dot product of
    Y(:,I) and Y(:,J) where Y(K,I) = -log(mean(X(K,:)) if X(K,I) = 1 and
    = log(mean(X(K,:)) if X(K,I) = 0. 
3. combineKernels.m
    given a MATLAB structure with each entry as an association network and a gene list, this function 
    extends the association networks so that they all include all the genes in the gene list. After calling
    this function, the genes in all association networks are ordered according to gene list. 
4. normalizeKernel.m
    Given an association network (kernel) it returs the normalize association network