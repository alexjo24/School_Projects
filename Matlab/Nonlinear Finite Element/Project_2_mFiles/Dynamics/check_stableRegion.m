clc
clear all
close all


% Gamma variation
gamma = [0:0.1:3];

%y=kx+m function
beta = (1/2)*gamma;

%Plot gamma vs beta.
figure
plot(gamma,beta,'LineWidth',1.5)
hold on
plot([1/2,1/2],[0 2],'LineWidth',1.5)
plot([0 3],[0.25,0.25], 'b--')
title('Stability')
xlabel('\gamma','FontSize', 14)
ylabel('\beta','FontSize', 14)
set(gca,'FontSize',14)
grid on