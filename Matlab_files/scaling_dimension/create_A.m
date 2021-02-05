function [ A ] = create_A( n )
    A = zeros(n,n);
    A(1,1:2) = [2, -1];
    A(n,(n-1):n) = [-1, 2];
    for k = 2:(n-1)
        A(k,(k-1):(k+1))=[-1, 2, -1];
    end
    A = ((n+1)^2/(1^2))*A;
end

