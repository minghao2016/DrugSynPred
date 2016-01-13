X = cell2mat(Mono.H);
X = ((X - min(X(:))) ./ (max(X(:) - min(X(:)))));

[ii, jj] = find(~isnan(X));
vv = X(sub2ind(size(X), ii, jj));
ksdensity(vv)

plot(sort(vv))


k = 10;
length = round(numel(vv) / k);
clear acc
for i = 1:k
    fprintf('%d\n', i)
    Y = X;
    mask_ii = ii((i-1)*length+1:i*length);
    mask_jj = jj((i-1)*length+1:i*length);
    mask_vv = vv((i-1)*length+1:i*length);
    Y(sub2ind(size(Y), mask_ii, mask_jj)) = nan;
    Z = cell2mat(fillinBlanks(num2cell(Y), C2C, D2D));
    
    acc(i) =  100*nanmean(abs(Z(sub2ind(size(Y), mask_ii, mask_jj)) - mask_vv));
    
end

mean(acc)