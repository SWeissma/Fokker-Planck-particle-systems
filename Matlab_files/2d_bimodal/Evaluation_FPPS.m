%%% Evaluation of Fokker--Planck Particle System

clear;
% load the corresponding data
% note that one has to run FPPS.m first
load('reference_data.mat')
load('ode_FPPS.mat')

% computation of the potentials
V_t = zeros(1,length(t));
for tk = 1:length(t)
    Xhelp = reshape(X(tk,:),[N_x,M]);
    V = zeros(M,1);
    for l = 1:M
        k = mvnpdf(Xhelp',Xhelp(:,l)',B)';
        V(l,:) = log(mean(k))+1/2*norm(sqrtm(R)\(h(Xhelp(:,l))-y),2)^2+1/2*norm(sqrtm(P0)\(Xhelp(:,l)-m0),2)^2;
    end
    V_t(1,tk) = mean(V);
end

V_t2 = zeros(1,length(t2));
for tk = 1:length(t2)
    Xhelp = reshape(X2(tk,:),[N_x,M]);
    V = zeros(M,1);
    for l = 1:M
        k = mvnpdf(Xhelp',Xhelp(:,l)',B)';
        V(l,:) = log(mean(k))+1/2*norm(sqrtm(R)\(h(Xhelp(:,l))-y),2)^2+1/2*norm(sqrtm(P0)\(Xhelp(:,l)-m0),2)^2;
    end
    V_t2(1,tk) = mean(V);
end

V_t3 = zeros(1,length(t3));
for tk = 1:length(t3)
    Xhelp = reshape(X3(tk,:),[N_x,M]);
    V = zeros(M,1);
    for l = 1:M
        k = mvnpdf(Xhelp',Xhelp(:,l)',B)';
        V(l,:) = log(mean(k))+1/2*norm(sqrtm(R)\(h(Xhelp(:,l))-y),2)^2+1/2*norm(sqrtm(P0)\(Xhelp(:,l)-m0),2)^2;
    end
    V_t3(1,tk) = mean(V);
end

% plot the potential
fig1 = figure(1);
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);

plot(t2,V_t2,'-.','Color',[0 0 0.7],'LineWidth',2,'DisplayName','preconditioned');hold on
plot(t,V_t,'-','Color',[0.7 0 0],'LineWidth',2,'DisplayName','derivative-free');hold on
plot(t3,V_t3,':','Color',[0 0.7 0],'LineWidth',2,'DisplayName','derivative-free localised');hold on
str = sprintf('potential');
xlabel('time','FontSize',20)
title(str,'FontSize',24)
legend('show','Location','northeast','FontSize',16,'Interpreter','latex')



% compute the kernel density estimate, weight function and variational
% derivative

Xhelp = reshape(X(end,:),[N_x,M]);
Xhelp2 = reshape(X2(end,:),[N_x,M]);
Xhelp3 = reshape(X3(end,:),[N_x,M]);

k = zeros(M,M);
for l = 1:M
    k(l,:) = mvnpdf(Xhelp',Xhelp(:,l)',B)';
end

h2 = zeros(M,M);
for l = 1:M
    h2(l,:) = mvnpdf(Xhelp2',Xhelp2(:,l)',B)';
end

h3 = zeros(M,M);
for l = 1:M
    h3(l,:) = mvnpdf(Xhelp3',Xhelp3(:,l)',B)';
end

F = zeros(100,100);
F2 = zeros(100,100);
F3 = zeros(100,100);
D = zeros(100,100);
D2 = zeros(100,100);
D3 = zeros(100,100);
Weights = zeros(100,100);
Weights2 = zeros(100,100);
Weights3 = zeros(100,100);

xtest = linspace(-3,3,100);
for k1 = 1:100
    for k2 = 1:100
        hx = mvnpdf(Xhelp',[xtest(k1),xtest(k2)],B)';
        hx2 = mvnpdf(Xhelp2',[xtest(k1),xtest(k2)],B)';
        hx3 = mvnpdf(Xhelp3',[xtest(k1),xtest(k2)],B)';
        F(k1,k2) = log(mean(hx))+1/2*norm(R^(1/2)\(h([xtest(k1);xtest(k2)])-y),2)^2+1/2*norm(P0\([xtest(k1);xtest(k2)]-m0),2)^2+sum(hx./sum(k));
        F2(k1,k2) = log(mean(hx2))+1/2*norm(R^(1/2)\(h([xtest(k1);xtest(k2)])-y),2)^2+1/2*norm(P0\([xtest(k1);xtest(k2)]-m0),2)^2+sum(hx2./sum(h2));
        F3(k1,k2) = log(mean(hx3))+1/2*norm(R^(1/2)\(h([xtest(k1);xtest(k2)])-y),2)^2+1/2*norm(P0\([xtest(k1);xtest(k2)]-m0),2)^2+sum(hx3./sum(h3));
        D(k1,k2) = mean(hx)*exp(sum(hx./sum(k)));
        D2(k1,k2) = mean(hx2)*exp(sum(hx2./sum(h2)));
        D3(k1,k2) = mean(hx3)*exp(sum(hx3./sum(h3)));
        Weights(k1,k2) = exp(sum(hx./sum(k)));
        Weights2(k1,k2) = exp(sum(hx2./sum(h2)));
        Weights3(k1,k2) = exp(sum(hx3./sum(h3)));
    end
end


% plot the variational derivative and the weight function
fig2 = figure(2);
clf(fig2)
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.9]);

[X,Y] = meshgrid(xtest,xtest);
subplot(3,2,1)
contourf(X,Y,F2,'edgecolor','none','DisplayName','preconditioned');
legend('show','Location','northeast','FontSize',16)
colorbar
str = sprintf('(a) variational derivative');
title(str,'FontSize',20)

subplot(3,2,2)
contourf(X,Y,Weights2,'edgecolor','none','DisplayName','preconditioned');
legend('show','Location','northeast','FontSize',16)
colorbar
str = sprintf('(b) implied weights');
title(str,'FontSize',20)

subplot(3,2,3)
contourf(X,Y,F,'edgecolor','none','DisplayName','derivative-free');
legend('show','Location','northeast','FontSize',16)
colorbar
str = sprintf('(c) variational derivative');
title(str,'FontSize',20)

subplot(3,2,4)
contourf(X,Y,Weights,'edgecolor','none','DisplayName','derivative-free');
legend('show','Location','northeast','FontSize',16)
colorbar
str = sprintf('(d) implied weights');
title(str,'FontSize',20)

subplot(3,2,5)
contourf(X,Y,F3,'edgecolor','none','DisplayName','derivative-free localised');
legend('show','Location','northeast','FontSize',16)
colorbar
str = sprintf('(e) variational derivative');
title(str,'FontSize',20)

subplot(3,2,6)
contourf(X,Y,Weights3,'edgecolor','none','DisplayName','derivative-free localised');
legend('show','Location','northeast','FontSize',16)
colorbar
str = sprintf('(f) implied weights');
title(str,'FontSize',20)


% plot of the kernel density estimate
fig3 = figure(3);
clf(fig3)
set(fig3, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.6]);
[X,Y] = meshgrid(xtest,xtest);
subplot(2,2,1)
contourf(X,Y,Cont','edgecolor','none','DisplayName','exact');

str = sprintf('(a) posterior density');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)

subplot(2,2,2)
contourf(X,Y,D2','edgecolor','none','DisplayName','preconditioned');
str = sprintf('(b) kernel density estimate');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)

subplot(2,2,3)
contourf(X,Y,D','edgecolor','none','DisplayName','derivative-free');
str = sprintf('(c) kernel density estimate');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)

subplot(2,2,4)
contourf(X,Y,D3','edgecolor','none','DisplayName','derivative-free localised');
str = sprintf('(d) kernel density estimate');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)
