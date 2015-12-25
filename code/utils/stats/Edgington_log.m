function [ combined_p ] = Edgington_log( p )
% Edgington method of combing p-values when the sum of p-values is less than or equal to 1: Edgington (1972) An additive method for combining probability values from independent experiments
    S = sum(p);
    K = numel(p);
    if(S <= 1)
        combined_p = S^K / factorial(K);
    else
        entries = arrayfun(@(j) (-1)^j * (S-j)^K / (factorial(K-j)*factorial(j)), 0:floor(S));
        combined_p = sum(entries);        
    end    
    if(combined_p <0 || combined_p > 1)
        combined_p = 1;
    end
end

