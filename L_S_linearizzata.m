%% LOOP SHAPING e SINTESI REGOLATORE
clear all;
close all;
clc;

T_a_5 = 0.03;
S = 0.02;
mu = 399;

xi = (abs(log(S)))/(sqrt(pi^2+log(S)^2))+0.2;
omega_c = 3/(xi*T_a_5);

s = tf('s');
F = 1/(1+2*xi*s/omega_c + s^2/omega_c^2);
L=mu*F/(1-F);
figure
hold on;
omega_min = 1;
omega_plot_min = 1e-1;
omega_plot_max = 1e3;
x_wc = [omega_plot_min; omega_plot_min; omega_c; omega_c];
y_wc = [-60; 0; 0; -60];
x_d = [omega_plot_min; omega_plot_min; 0.5; 0.5];
y_d = [40; 0; 0; 40];
x_n = [1e5; 1e5; 5e6; 5e6];
y_n = [-63; 200; 200; -63];
patch(x_wc, y_wc, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_d, y_d, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_n, y_n, 'r', 'FaceAlpha', 0.1);
bode(L);

figure;
YY_w = step(F);
%YY_d = step(1/(1+L));
%YY = YY_w + YY_d;
plot(YY_w);