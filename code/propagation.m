% Propagate initial vector through network
% network is already column_normalized
function prop_vector = propagation(network , alpha, initial)

THRESHOLD = 1e-10;
MAXIT = 100;
n = size(network, 1);
L = network - diag(sum(network,2));
L = eye(n) - alpha*L;
coef = inv(L);

residue = 1;
iter = 1;

prop_vector = initial;

while (residue > THRESHOLD && iter < MAXIT),
    old_prop_vector = prop_vector;
    prop_vector     = coef*prop_vector;
%     prox_vector     = prox_vector/norm(prox_vector, 1);
    residue         = norm(prop_vector-old_prop_vector);
    iter            = iter + 1; 
end


end
