clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LIMITI PLOT
omega_plot_min = 1e-2;
omega_plot_max = 5*1e6;
degree_plot_min = -360;
degree_plot_max = 360;
dB_plot_min = -2*1e2;
db_plot_max = 2*1e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETRI DATI
k = 100; % elasticità del disco
beta = 0.6; % attrito viscoso
alpha = 30; % angolo in gradi tra i due alberi
J = 800; % momento di inerzia della tavola
theta_e = 120; % posizione angolare della tavola di equilibrio

alpha_rad = deg2rad(alpha);
theta_e_rad = deg2rad(theta_e);

tau = @(theta) cos(alpha_rad)/(1-(sin(alpha_rad)*cos(theta))^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% COPPIA DI EQUILIBRIO E MATRICI DEL SISTEMA LINEARE (PUNTO 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% COPPIA DI EQUILIBRIO
x_e1 = theta_e_rad;
x_e2 = 0;
u_e  = k*x_e1/tau(x_e1);
x_e  = [x_e1;x_e2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MATRICI A,B,C,D
df1_dx1_e = 0;
df1_dx2_e = 1;
df2_dx1_e = -2*u_e*(cos(alpha_rad) * sin(alpha_rad)^2 * sin(x_e1) * cos(x_e1))/(J * (1 - sin(alpha_rad)^2 * cos(x_e1)^2)^2 )-k/J;
df2_dx2_e = -beta/J;

A = [df1_dx1_e df1_dx2_e; df2_dx1_e df2_dx2_e];
B = [0; tau(x_e1)/J];
C = [1 0];
D = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINIZIONE DELLA G (PUNTO 2)
s  = tf('s');
GG = C*inv(s*eye(2) - A)*B + D; %GG = 0.001155/(s^2+0.00075*s+0.06454);
margin(GG,{omega_plot_min,omega_plot_max});

legend_arg = ["G(j\omega)"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ZONE PROIBITE (PUNTO 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPECIFICHE DATE
%% e_inf <= 0.005, w(t) = 1(t) e d(t) = 1(t)
%% M_f >= 30 gradi
%% S% <= 2%
%% Ta,5% < 0.03s
%% d(t) con banda [0, 0.5], A_d = 40dB
%% n(t) con banda [1e5, 5*1e6], A_n = 63dB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINIZIONE SPECIFICHE
e_star = 0.005;
WW = 1;
DD = 1;
Mf_min = 30;
S_star = 0.02;
T_a_5_star = 0.03;
A_d = 40;
omega_d_min = 0;
omega_d_max = 0.5;
A_n = 63;
omega_n_min = 1e5;
omega_n_max = 5*1e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CALCOLO PARAMETRI PER PATCHES
xi_min = abs(log(S_star))/sqrt(pi^2 + log(S_star)^2); %+ 0.01
Mf_star = max(100*xi_min, Mf_min);
omega_c_min = 3 * 100 / (Mf_star * T_a_5_star);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VERTICI PATCH DISTURBO SULL' USCITA
x_d = [omega_plot_min; omega_plot_min; omega_d_max; omega_d_max];
y_d = [dB_plot_min; A_d; A_d; dB_plot_min];

%% VERTICI PATCH DISTURBO DI MISURA
x_n = [omega_n_min; omega_n_min; omega_n_max; omega_n_max];
y_n = [-A_n; omega_plot_max; omega_plot_max; -A_n];

%% VERTICI PATCH PULSAZIONE CRITICA
x_c = [omega_plot_min; omega_plot_min; omega_c_min; omega_c_min];
y_c = [dB_plot_min; 0; 0; dB_plot_min];

%% VERTICI PATCH MARGINE DI FASE
x_f = [omega_c_min; omega_c_min; omega_n_min; omega_n_min];
y_f = [degree_plot_min; Mf_star-180; Mf_star-180; degree_plot_min];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PATCH E PLOT DI G DOMINIO FREQUENZE
figure;
patch(x_d, y_d, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_n, y_n, 'y', 'FaceAlpha', 0.1);
hold on;
patch(x_c, y_c, 'b', 'FaceAlpha', 0.1);
hold on;
margin(GG,{omega_plot_min,omega_plot_max});
grid on; zoom on;
hold on;
patch(x_f,y_f, 'g', 'FaceAlpha', 0.1);

legend_arg = ["\omega_{c_{min}}", "A_n", "A_d", "G(j\omega)", "Mf"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PATCH E PLOT G DOMINIO TEMPO
% FF=GG/(1+GG);
% T_simulation = 2*T_a_5_star;
% LV = evalfr(WW*FF,0);
% 
% [y_w,t_step] = step(FF,T_simulation);
% 
% figure;
% plot(t_step,y_w,'g');
% hold on;grid on;
% 
% % vincolo sovraelongazione
% patch([0,T_simulation,T_simulation,0],[LV*(1+S_star),LV*(1+S_star),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
% 
% % vincolo tempo di assestamento all'5%
% patch([T_a_5_star,T_simulation,T_simulation,T_a_5_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
% patch([T_a_5_star,T_simulation,T_simulation,T_a_5_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
% 
% legend_arg = ["y_w(t)", "S%", "T_{a,5%}"];
% legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOOP SHAPING E SINTESI REGOLATORE (PUNTO 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SINTESI REGOLATORE STATICO
G_0 = abs(evalfr(GG, j * 0));
mu_min_e = (WW + DD)/e_star;
mu_min_d = 10^(A_d/20);
mu_s_min = max(mu_min_e, mu_min_d)/G_0;
mu_s = mu_s_min ;
RR_s = mu_s; %regolatore statico
GG_e = RR_s * GG;

figure;
patch(x_d, y_d, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_n, y_n, 'y', 'FaceAlpha', 0.1);
hold on;
patch(x_c, y_c, 'b', 'FaceAlpha', 0.1);
hold on;
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;
hold on;
patch(x_f,y_f, 'g', 'FaceAlpha', 0.1);

legend_arg = ["\omega_{c_{min}}", "A_n", "A_d", "G_e(j\omega)", "Mf"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SINTESI REGOLATORE DINAMICO

omega_c_star = 200; % scegliamo omega_c_star

mag_omega_c_star = abs(evalfr(GG_e,j*omega_c_star));
arg_omega_c_star = rad2deg(angle(evalfr(GG_e,j*omega_c_star)));

M_star = 1/mag_omega_c_star;
phi_star = Mf_star + 5 - 180 - arg_omega_c_star;

tau = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha = alpha_tau / tau;

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi');
    return;
end
%R_high_frequency = 1/(1 + s/5e4); % per robustezza nell'attenuazione del disturbo di misura preferiamo inserire un polo ad alte frequenze e abbattere maggiormente
RR_d = (1 + tau*s)/(1 + alpha * tau*s);%*R_high_frequency; % regolatore dinamico
RR = RR_s*RR_d;
LL = RR_d * GG_e; % funzione d'anello finale

figure;
patch(x_d, y_d, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_n, y_n, 'y', 'FaceAlpha', 0.1);
hold on;
patch(x_c, y_c, 'b', 'FaceAlpha', 0.1);
hold on;
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;
hold on;
patch(x_f,y_f, 'g', 'FaceAlpha', 0.1);

legend_arg = ["\omega_{c_{min}}", "A_n", "A_d", "L(j\omega)", "Mf"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TEST DEL SISTEMA LINEARIZZATO NEL TEMPO (PUNTO 4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINIZIONE DELLE FUNZIONI DI SENSIBILITA'
FF = LL/(1+LL);
SS = 1/(1+LL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RISPOSTA A w(t)=1(t)

T_simulation = 2*T_a_5_star;
LV = evalfr(WW*FF,0);

[y_w,t_step] = step(FF,T_simulation);

figure;
plot(t_step,y_w,'g');
hold on;grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PATCHES ZONE PROIBITE DOMINIO DEL TEMPO

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star),LV*(1+S_star),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
% vincolo tempo di assestamento all'5%
patch([T_a_5_star,T_simulation,T_simulation,T_a_5_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a_5_star,T_simulation,T_simulation,T_a_5_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

legend_arg = ["y_w(t)", "S%", "T_{a,5%}"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RISPOSTA A d(t)='serie di fourier'

figure;
dd = 0;
tt = 0:1e-2:2e2;
for k=1:4
    dd = 0.2*sin(0.1*k*tt) + dd;
end
y_d = lsim(SS,dd,tt);
plot(tt,dd,'r')
hold on;
plot(tt,y_d,'g');
grid on;

legend_arg = ["d(t)", "y_d(t)"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RISPOSTA A n(t)='serie di fourier'

figure;
nn = 0;
tt = 0:1e-5:2*1e-3;
for k=1:4
    nn = 0.2*sin(1e5*k*tt) + nn;
end
y_n = lsim(-FF,nn,tt);
plot(tt,nn,'r')
hold on;
plot(tt,y_n,'g')
grid on;

legend_arg = ["n(t)", "y_n(t)"];
legend(legend_arg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%