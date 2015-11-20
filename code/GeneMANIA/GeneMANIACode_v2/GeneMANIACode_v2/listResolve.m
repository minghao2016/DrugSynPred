%LISTRESOLVE - finds the indices of a sublist in a reference list
%
%  IX = listResolve(REF_LIST, SUB_LIST)
%
%  IX(I) is the index of SUB_LIST{I} in REF_LIST.
%
%  If SUB_LIST(I) does not correspond to anything in REF_LIST,
%  or if SUB_LIST(I) appears previously in SUB_LIST, IX(I) = NaN

function ix = listResolve(reflist, sublist)

[overlap, ia, ib] = intersect(reflist, sublist);
%if length(overlap) ~= length(sublist)
%  warning('Second argument is not a sublist of the first');
%end

ia(end+1) = NaN;
[junk, order] = sort(ib);
map = repmat(length(ia), 1, length(sublist));
map(junk) = order;
ix = ia(map);
