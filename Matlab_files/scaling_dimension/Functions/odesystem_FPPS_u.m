function [ f ] = odesystem_FPPS_u(X,u,N_x,R,y,m0,P0,delta)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Input:
    % X = current state of the ODE
    % u = forward model (KL expansion)
    % N_x = dimension of the state
    % R = noise covariance
    % y = observation
    % m0 = prior mean
    % P0 = prior covariance
    % delta = kernel scaling parameter
    %%% Output:
    % f = drift term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = length(X)/(N_x);
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % computation of the sample covariances
    xquer = mean(Xhelp,2);
    uquer = mean(u(Xhelp),2);
    P_t = 1/M*(Xhelp-xquer)*(Xhelp-xquer)';
    Pu_t = 1/M*(Xhelp-xquer)*(u(Xhelp)-uquer)';
    
    B_est = diag(P_t);
    
    % optimal kernel bandwidth wrt. AMISE
    kernel_width = ((4/((2+N_x)*M))^delta);
    
    % adapting the kernel covariance
    B = kernel_width*diag(B_est);

    % computation of the kernels 
    k = zeros(M,M);
    for j = 1:M
        k(j,:) = (1/sqrt(2*pi))^N_x*1/sqrt(det(B))*exp(-1/2*sum((sqrtm(B)\(Xhelp-Xhelp(:,j))).^2,1));
    end    
    
    %computation of the drift
    K1 = zeros(N_x,M);
    K2 = zeros(N_x,M);
    for i=1:M
        N1 = (B\(Xhelp(:,i)-Xhelp)).*k(i,:);
        K1(:,i) = sum(N1,2)/sum(k(i,:));
        K2(:,i) = sum(N1./sum(k,1),2);
%         K2(:,i) = mean(Norm1.*h_j,2);
    end  
    f = -Pu_t*((R\(u(Xhelp)-y)))-P_t*(P0\(Xhelp-m0))+P_t*((K1)+(K2));  
    f = f(:);
    
end