function Linv = getInverseLaplacian(Wk)
% function Linv = getInverseLaplacian(Wk)
%
%   Wk -- an affinity network represented as a kernel
% Compute D^(-1/2) L D^(-1/2), L = D - Wk, turns into this with some algebra.


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
Linv = inv(M);
