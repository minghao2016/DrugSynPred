function ix = lookup(list, ref)
%LOOKUP - lookup a list in another one
%
%  INDICES = LOOKUP(LIST, REF)
%
%   Looks up LIST in REF, INDICES(I) is the index of LIST(I) in
%   REF, INDICES(I) = NaN if LIST is not in REF.  

[uniq, junk, ixs] = unique(list);
map = listResolve(ref, uniq);
ix = map(ixs);



