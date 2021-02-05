function [ f ] = odesystem_FPPS_fix(X, h, y, R, N_x,P0,m0,B)
    
    M = length(X)/(N_x);
    
    % collecting the current ensemble of particles
    Xhelp = reshape(X,[N_x,M]);
    
    % computation of the sample covariances
    hx = h(Xhelp);
    huquer = mean(hx,2);
    uquer = mean(Xhelp,2);
    Pxx = 1/M*Xhelp*Xhelp'-uquer*uquer';
    Pxh = 1/M*Xhelp*hx'-uquer*huquer';
    

    % computation of the kernels 
    k = zeros(M,M);
    for j = 1:M
        k(j,:) = (1/sqrt(2*pi))^N_x*1/sqrt(det(B))*exp(-1/2*sum((sqrtm(B)\(Xhelp-Xhelp(:,j))).^2,1));
    end
    
    % computation of the drift
    K1 = zeros(N_x,M);
    K2 = zeros(N_x,M);
    for i=1:M
        N1 = (Xhelp(:,i)-Xhelp).*k(i,:);
        K1(:,i) = sum(N1,2)/sum(k(i,:));
        K2(:,i) = sum(N1./sum(k,1),2);
    end  
    f = (Pxh*(R\(y-hx))-Pxx*(P0\(Xhelp-m0)))+Pxx*(B\K1)+Pxx*(B\K2);    
    f = f(:);
    
    
end

