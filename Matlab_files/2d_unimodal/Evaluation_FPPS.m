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

% plot the potential
fig1 = figure(1);
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);

plot(t,V_t,'-','Color',[0.7 0 0],'LineWidth',2,'DisplayName','derivative-free');hold on
str = sprintf('potential');
xlabel('time','FontSize',20)
title(str,'FontSize',24)
legend('show','Location','northeast','FontSize',16,'Interpreter','latex')



% compute the kernel density estimate

Xhelp = reshape(X(end,:),[N_x,M]);

k = zeros(M,M);
for l = 1:M
    k(l,:) = mvnpdf(Xhelp',Xhelp(:,l)',B)';
end

D = zeros(100,100);

xtest = linspace(-2,2,100);
for k1 = 1:100
    for k2 = 1:100
        hx = mvnpdf(Xhelp',[xtest(k1),xtest(k2)],B)';
        D(k1,k2) = mean(hx)*exp(sum(hx./sum(k)));
    end
end


% plot of the kernel density estimate
fig2 = figure(2);
clf(fig2)
set(fig2, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);
[X,Y] = meshgrid(xtest,xtest);
subplot(1,2,1)
contourf(X,Y,Cont','edgecolor','none','DisplayName','exact');

str = sprintf('(a) posterior density');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)

subplot(1,2,2)
contourf(X,Y,D','edgecolor','none','DisplayName','derivative-free');
str = sprintf('(b) kernel density estimate');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)
