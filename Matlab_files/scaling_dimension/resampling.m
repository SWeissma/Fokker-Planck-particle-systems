function [Z_new] = resampling(Z,X,B,N_new)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % Z = sample
    % X = particles/ nodes for the kernel
    % B = Gaussian kernel covariance
    % N_new = new number of samples
    %%% Output:
    % Z_new = resampled sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [N_x,N] = size(Z);
    w = zeros(1,N);
    [~,M] = size(X);
    
    % computation of the weights
    k = zeros(M,M);
    for l = 1:M
        k(l,:) = 1/sqrt(det(B))*exp(-1/2*sum((sqrtm(B)'\(X-X(:,l))).^2,1));
    end
    
    for l = 1:N
        kx = 1/sqrt(det(B))*exp(-1/2*sum((sqrtm(B)'\(X-Z(:,l))).^2,1));
        w(l) = exp(mean(kx./mean(k,1)));
    end
    
    w = w./sum(w);
    
    % resampling
    Z_new = zeros(N_x,N_new);
    
    for n = 1:N_new
        X = sum(w)*rand;
        Q = cumsum(w)-X;
        Q = Q(Q<0);
        i = length(Q)+1;
        Z_new(:,n) = Z(:,i);
    end
end

