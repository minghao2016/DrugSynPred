function [f,numIt] = getScoreVectorCG(Y, Wk)
%getScoreVector -- Compute score vector given a kernel and a label vector
%                  using the conjugate gradient method
%
%   F = getScoreVector(Y, Wk) 
%
%   Wk -- an affinity network represented as a kernel
%   Y  -- label vector: -1/+1 for training points, 0 for other points
%
% Returns F, a vector of discriminant values ("scores") derived by solving
% the system Z = Wk * F. where Z(I) = Y(I) if Y(I) = -1/+1, otherwise
%  Z(I) = (N-M)/(N+M) where N = nnz(Y = 1) and M = nnz(Y = -1), numIt is
%  the number of conjugate gradient iteration until convergence
%
% $Revision: 1.4 $

% Compute D^(-1/2), applying a tiny bit of regularization to make sure
% we're invertible.

% Compute D^(-1/2) L D^(-1/2), L = D - Wk, turns into this with some algebra.


Y = Y(:);
NN = nnz(Y == 1);
MM = nnz(Y == -1);
Y(find(Y==0)) = (NN-MM) / (NN + MM); % set the label bias to the mean of the labels

twoI = sparse(1:length(Wk), 1:length(Wk), 2);
oneI = sparse(1:length(Wk), 1:length(Wk), 1);

d = 1 ./ sqrt(sum(Wk,1) + eps);
[rr,cc,ss] = find(Wk);
Shat = ss(:)' .* d(rr) .* d(cc);
lD = length(d);
What =  sparse(rr,cc,Shat,lD,lD);
Dhat = sparse(1:lD, 1:lD, sum(What, 1));
Lhat = Dhat - What;
M = oneI + Lhat;
M = max(M,M'); % make sure M is symmetric

finit = zeros(size(Y,1),1);
maxit = 1000;  % hard-coded maximum number of iterations to 1000: usually it takes less than 25 iterations
[f,res,numIt] = conjGrad(finit,Y,M,maxit);

