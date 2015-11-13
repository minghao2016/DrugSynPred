function WWW = normalizeKernel(WW)
% function WWW = normalizeKernel(WW)

% Subtract any diagonal scores
WW = WW - diag(diag(WW));

% Get the indices and values of all the non-zero elements
[rr, cc, ss] = find(WW);
% Get the column sums, which are the same as the row sums since the matrix is symmetrical
DD = full(sum(WW, 1));
% Create a filehandle to represent the normalization equation
norm = @(x) (1./sqrt(x));
% Only normalize for the rows/columns where there is at least one non-zero interaction
normalize = spfun(norm, DD);
% Normalize by dividing the value by the product of the square-root of the column sum and the square-root of the row sum
WWW = sparse(rr, cc, ss' .* normalize(rr) .* normalize(cc), size(WW, 1), size(WW,2));