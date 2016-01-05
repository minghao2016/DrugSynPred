% running TIMMA test for different cell lines.
% this is the function to find the optimal find by k_max, S ,y

function[minInd] = TIMMATest2(S, y, k_max)
[d,k] = size(S);
error = zeros(1,k);
for i= 1:23
    disp(i)
    [M,err,kset] = TIMMA_floating2(S,y,i,k_max,1);
    error(i) = mean(err);
end
[~, minInd] = min(error);

test=S(:,7);
tmp2 = [];
for i=1:101
    tmp = sum(abs(test - S(:,i)));
    if tmp==0
        tmp2 = [tmp2 i];
    end
end