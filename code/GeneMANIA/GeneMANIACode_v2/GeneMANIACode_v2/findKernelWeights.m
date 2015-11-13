function [alpha,indices] = findKernelWeights(y, W, myKtK)
%findKernelWeights -- find least squares estimates for kernel weights
%
%  [ALPHA, INDICES] = findKernelWeights(y, W)
%
%   y -- vector of labels  y(i) is -1, 0, or 1
%   W -- cell array of sparse matrices.  W{I} is the I-th affinity network
%        aka weight matrix aka kernel.  W{I} is assumed to be
%        sparse and symmetric.
%  
%  If calculations yield negative weights for any kernels, those kernels
%  are removed and alpha is recomputed with the remaining kernels. This 
%  process is repeated until only positive weights remain.
%  This function eliminates terms which correspond to negative/negative 
%
%  pairs in the target matrix before performing the linear regression. 
%  Thus far it appears that this modification does not result in significantly
%  better results, but may speed up computation significantly in some
%  situations.
%
%  Returns ALPHA, a vector of weights for matrices (not including those 
%  which received 0 weight), and INDICES, a vector containing the indices 
%  of the kernels which had non-zero weights, in the same order as the
%  elements of ALPHA.
%
%  If and ONLY IF a third output argument is specified, the regression is
%  done with a bias term.
%  $ Revision 1.2 $

% ---- Setting a few constants
ixPos = find(y == 1);
ixNeg = find(y == -1);
nPos = nnz(y == 1);
nNeg = nnz(y == -1);
NN = length(y);
MM = length(W);
usesBias = 1 ; % always use bias
epsilon = eps;             % our numerical def'n of zero

posConst = 2*nNeg / (nPos + nNeg);   % Label bias for +/+ elements
negConst = -2*nPos / (nPos + nNeg);  % Label bias for +/- elements

% ---- A few checks on the arguments.

if nPos == 0 | nNeg == 0
    error('Must have at least one positive and one negative');
end

nZero = nnz(y == 0);
if nPos + nNeg + nZero ~= length(y)
    error('Elements of Y must be either -1, 0, or 1');
end

if NN == 0, error('Empty target vector!'); end
if MM == 0, error('Empty kernel set'); end

% ---- Calculation of the linear indices in W{K} of all the 
%      (i,j)-elements for which
%  ixpp:  y(i) = 1 and y(j) = 1 and i ~= j,   (aka +/+)
%  ixpn:  y(i) = 1 and y(j) = -1              (aka +/-)

% Note: because the affinity matrices in W are symmetric, all the sums
% and dot products involving +/- (i,j) elements are the same for the
% -/+ (j,i) elements.  So, we simply calcuate double any sums and dot
% products that we calculate over the +/- elements.

% ---- Calculation of the two matrices needed to do the linear regression:
%  KtK:  the gram matrix containing the dot products of the kernels
%  KtT:  the dot products of the kernels with the target values

KtK = zeros(MM+1, MM+1);
KtT = zeros(MM+1, 1);

Wpp = cell(1, MM);    % temp storage of +/+ non-diagonal affinities
Wpn = cell(1, MM);    % temp storage of +/- elements affinities
nPpElem = nPos * (nPos - 1); % # +/+ elements that aren't diagonal
nPnElem = 2 * nPos * nNeg;   % # +/- elements
ppTarget = posConst^2;           % target for +/+ elements
pnTarget = posConst * negConst;  % target for +/- elements
elem = (nPpElem + nPnElem);
% value of the bias interactions
biasVal = 1 / (nPpElem + nPnElem);

KtT(1) = biasVal * (ppTarget * nPpElem + pnTarget * nPnElem);
KtK(1,1) = biasVal;

for ii = 1:MM
  Wpp{ii} = W{ii}(ixPos, ixPos);
  tmp = sparse(1:nPos, 1:nPos, diag(Wpp{ii}));
  Wpp{ii} = Wpp{ii} - tmp;   % sets the diagonal to zero
  Wpn{ii} = W{ii}(ixPos, ixNeg);

  ssWpp = full(sum(sum(Wpp{ii})));
  ssWpn = full(sum(sum(Wpn{ii})));
  
  KtT(ii+1) = ppTarget * ssWpp + 2 * pnTarget * ssWpn;
  KtK(ii+1, 1) = biasVal * (ssWpp + 2 * ssWpn);
  KtK(1, ii+1) = KtK(ii+1, 1);

  for jj = 1:ii
    KtK(ii+1,jj+1) = full(...
        sum(sum(Wpp{ii} .* Wpp{jj})) + ...
        2 * sum(sum(Wpn{ii} .* Wpn{jj})) ...
        );
    KtK(jj+1,ii+1) = KtK(ii+1,jj+1);   % make KtK symmetric
  end
end

% ---- Removes the bias term if its not used

if ~usesBias
  KtT = KtT(2:end);
  KtK = KtK(2:end, 2:end);
end

% ---- Removing empty columns

ss = sum(abs(KtK), 2);
goodIdx = find(sum(abs(KtK), 2) > epsilon * max(ss));
KtK = KtK(goodIdx, goodIdx);
KtT = KtT(goodIdx);

indices = 1:(length(goodIdx) - 1);

done = 0;


while ~done,

    % alpha = (K' * K)^-1 K' * T(:)
    alpha = inv(KtK) * KtT;
    if(any(isnan(alpha)))
        indices = [];
        break
    end

    % Find the locations of all the negative weights
    negWeights = find(alpha < 0);

    negWeights = negWeights(:);
    negWeights = setdiff(negWeights, 1);  % always have a bias term now
    if length(negWeights) > 0,

        % Save ourselves some computation by just eliminating rows and cols
        % from K^T K and elements from K^T y. The savings from this appear
        % to be nominal, though.
        KtK(:,negWeights) = [];
        KtK(negWeights,:) = [];
        KtT(negWeights) = [];

        % Remove from the indices array.
        indices((negWeights - 1)) = [];
    else,
        done = 1;
        fprintf('%d matrices chosen, indices %s, weights %s\n',...
            size(KtK,1)-1,mat2str(indices),mat2str(alpha(2:end)));
    end;
end;




% Undo our previous trickery
indices = goodIdx(indices);


% if we have no non-zero alpha entries, use the average kernel
if length(indices) == 0,
  indices = 1:length(W);
  alpha = [0 1/length(W)*ones(1,length(W))]';
  fprintf('All kernels eliminated or empty, using average kernel\n');
end;
% ignore the bias term which is alpha(1)
alpha = alpha(2:end);
