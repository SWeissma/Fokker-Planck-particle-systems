function [Z] = sampling_kernel(X,C,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = particles/ nodes for the kernel
    % C = Gaussian kernel covariance
    % N = number of samples
    %%% Output:
    % Z = generated sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [k1,M] = size(X);
    % sampling of each kernel
    Z = X(:,1)+mvnrnd(zeros(k1,1),C,N)';
    for m=2:M
        Z = [Z,X(:,m)+mvnrnd(zeros(k1,1),C,N)'];
    end
    %%% alternative sampling on sum over all kernels
    %     m = ceil(M*rand(1,N));
    %     Z = X(:,m)+mvnrnd(zeros(k1,1),C,N)';
end

