clc
clear all
close all

%Material parameters
E = 1;
A0 = 1;
v = 0.3;

G = E/(2*(1+v));

%Stretch , Lambda
Lambda = [0:0.01:1.4];

%Nonlinear material model
S = G*(Lambda - (1./(Lambda.^2)));
N = S*A0;
F = N.*Lambda;


%Green's strain
Lambda2 = [-1.4:0.01:1.4];
eg = (1/2).*((Lambda2.^2)-1);
S = E*eg;
N = S*A0;
F2 = N.*Lambda2;


figure
plot(Lambda, F,'LineWidth',1.5)
title('Force-Stretch in a bar element')
xlabel('Stretch')
ylabel('F [N]')
grid on
hold on
plot(Lambda2, F2,'LineWidth',1.5)
axis([-1.4 1.4 -3 3])
legend('True strain','Green´s strain')