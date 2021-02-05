%%% Evaluation Interacting Langevin particle dynamics

clear;
% load the corresponding data
% note that one has to run Langevin.m first
load('reference_data.mat')
load('SDE_Langevin.mat')

fig1 = figure(1);
clf(fig1)
set(fig1, 'Units', 'normalized', 'Position', [0.1, 0.5, 0.6, 0.3]);

subplot(1,2,1)
contourf(z_axis,z_axis,Cont','edgecolor','none','DisplayName','posterior distribution');hold on
str = sprintf('(a) exact');
title(str,'FontSize',20)
legend('show','Location','northeast','FontSize',16)

subplot(1,2,2)
Xest2 = reshape(Xest2_path(end,:),[N_x,M2]);
patch([-3 -3 3 3],[-3 3 3 -3],[61 38 168]./255);hold on
p = scatter(Xest2(1,:),Xest2(2,:),'MarkerFaceColor',[0 .9 0],'DisplayName','derivative-free');hold on
str = sprintf('(b) particle system');
title(str,'FontSize',20)
legend('show',p,'Location','northeast','FontSize',16)
