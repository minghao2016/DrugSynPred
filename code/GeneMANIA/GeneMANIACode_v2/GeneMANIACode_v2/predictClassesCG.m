function [p, k, wts] = predictClassesCG(labels, kernels, Regularize)
%predictClasses -- make predictions given a set of kernels and labels
%
%  [P, K, WTS] = predictClasses(LABELS, KERNELS)
%
%   LABELS    -- a vector of consisting of 1's or -1's for the labeled
%                elements, and 0's for the unlabeled elements
%   KERNELS   -- a cell array of kernel matrices (assumed to be sparse)
%   Regularize - optional, if specified (Regularize == 1), a regularized
%   version of GeneMANIA network integration is used. 
%   
%  The kernels are combined as a weighted sum, with weights determined by 
%  the MANIA algorithm, and a score vector is computed from this combined
%  kernel.
%
%  Returns P, a real-valued score vector that can be thresholded appropriately
%  (no thresholding is applied by this function), K,  a list of indices of 
%  the kernels which received nonzero weight, and WTS, the corresponding 
%  weights assigned to the kernels indexed in K (in the same order as K).
%
% by David Warde-Farley
% $Revision 1.1$ Sara: using the bias term is now mandatory
% $Revision 1.2$ Sara: an addition input parameter which indicates equal
% weight prior regularization

if length(labels) ~= length(kernels{1}),
    error('kernels must be same dimension as label vector');
end


if length(kernels) == 1 % if there is only one kernel then dont get kernel weights
    wts = 1;
    k = 1;
    bias = 1;
elseif nargin > 2
    [wts, k] = findKernelWeightsReg(labels, kernels);
else
    [wts, k] = findKernelWeights(labels, kernels);
end

nNonzero = 0;

for i = 1:length(k)
    nNonzero = nNonzero + nnz(kernels{k(i)});
end

combined = sparse([],[],[],length(labels),length(labels),nNonzero);

for i = 1:length(k),
    combined = combined + wts(i) *(kernels{k(i)});
end;

wts = wts(:);

combined = sparse(combined);
[p,numIt] = getScoreVectorCG(labels, combined);
