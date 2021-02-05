function [O] = observation_matrix(K,I)
    O = zeros(K,I-1);
    l = I/(K+1);
    for k = 1:(K)
        O(k,k*l) = 1;
    end
end

