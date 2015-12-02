% Propagate initial vector through network
% network is already column_normalized
function prop_vector = propagation(network , alpha, initial)

    n = size(network, 1);
    L = network - diag(sum(network,2));
    L = eye(n) - alpha*L;
    coef = inv(L);
 
    prop_vector = coef * initial ; 


end
