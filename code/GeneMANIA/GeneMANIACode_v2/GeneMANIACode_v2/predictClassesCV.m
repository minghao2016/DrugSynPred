function [areas, r, b, weights] = predictClassesCV(labels, kernels, nFolds, regularize, pp)
%predictClassesCV -- Perform cross-validation on MANIA
%
%  [areas, r, b, weights] = predictClassesCV(labels, kernels, nFolds, pp)
%  
%  LABELS    -- a label vector, 1 for positive examples, -1 for negative
%               examples, for unlabelled examples (ignored by this function)
%  KERNELS   -- A cell array of kernel matrices, assumed to be sparse
%  NFOLDS    -- The number of folds of cross-validation to perform
%  regularize - A scalar which specify which version of findKernelWeights
%               to use, regularize=1 uses equal weight prior otherwise 
%               un-regularized network integration is used.   
%  PP        -- Optional vector of length length(find(labels)) containing a
%               random permutation of the indices 1:length(find(labels)). 
%               You must specify a weighting scheme in the previous parameter
%               in order to use this argument.
%
% Returns AREAS, a vector of length NFOLDS containing the ROC area on each step
% of cross-validation. Optionally returns R, a cell array of length NFOLDS, 
% containing the score vectors for each cross-validation step, and B, which
% returns a score vector of leave-out scores for each gene (i.e., the score
% that gene got when it was one of the left out genes).
%
% $Revision: 1.2 $ Sara: Using bias is now the default
% $Revision: 1.3 $ Sara: Now users must specify if they want regularized or
% unregularized network integration. 

labelIx = find(labels);
NN = length(labelIx);

if (nargin < 5) % if permutation of the indices for cross-validation is not given
    rand('state', sum(100*clock));
    pp = randperm(NN);
end

offset = 0;

for ii = 1:nFolds
    lastElem = min(NN, offset+floor(NN/nFolds));
    ix = labelIx(pp(offset+1:lastElem));
    offset = lastElem;
    myLabels = labels;
    myLabels(ix) = 0;
    
    % using the bias is now hard-coded
    if (regularize)
        [res,k,wts] = predictClassesCG(myLabels, kernels,regularize);
    else
        [res,k,wts] = predictClassesCG(myLabels, kernels);
    end
    
    areas(ii) = calcROCarea(res(ix), round((1+labels(ix))/2));
    
    if nargout > 1
        r{ii} = res;
    end

    if nargout > 2,
        b(ix) = res(ix);
    end

    if nargout > 3,
        weights(ii,:) = zeros(1,length(kernels));
        weights(ii,[(k(:) )']) = wts;     
    end;
    
end

fprintf('areas on %d-fold cross-validatoin:  %s\n',...
            nFolds,mat2str(areas));