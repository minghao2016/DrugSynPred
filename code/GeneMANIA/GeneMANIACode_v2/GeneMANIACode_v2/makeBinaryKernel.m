function G = makeBinaryKernel(X, K)
%function G = makeZhangKernel(X, K)
%
%  X is a MxN matrix of binary data
%  Columns of X are the datapoints.
%  Returns a NxN sparse matrix G
%  G(I,J) = 0 if J is not within the K closest neighbours
%  otherwise, G(I,J) = "Fisher kernel" of X(:,I) and X(:,J)
%   in this case, the Fisher kernel is the dot product of
%   Y(:,I) and Y(:,J) where Y(K,I) = -log(mean(X(K,:)) if X(K,I) = 1 and
%                                  = log(mean(X(K,:)) if X(K,I) = 0

ss = sum(X, 2);
goodIx = find(ss);
if length(goodIx) < length(ss)
    X = X(goodIx, :);
end

[M,N] = size(X);
mm = full(mean(X, 2));
f0 = log(1-mm);
f1 = -log(mm);
tt = f0 .* (f1 - f0);
uu = (f1 - f0).^2;

zz = f0' * f0;

[rr,cc] = find(X);
T = sparse(rr,cc,tt(rr), M, N);
hh = sum(T,1);

G = sparse([],[],[],N,N,N*(K+1));
UU = sparse(cc, rr, uu(rr), N, M);
fprintf('makeKernel called with %d datapoints.  One dot / 1000 datapoints\n', N);

for ii = 1:N
    vv = UU(ii,:) * X;
    value = zz + hh(ii) + hh + vv;
    [sortedValues, indSort] = sort(value,'descend');
    G(indSort(1:K+1), ii) = sortedValues(1:K+1)';
    if mod(ii, 1000) == 0
        fprintf('.');
        if mod(ii, 8000) == 0
            fprintf('\n');
        end
    end
end
fprintf('\n');

G = G - diag(diag(G));
G = max(G, G');
