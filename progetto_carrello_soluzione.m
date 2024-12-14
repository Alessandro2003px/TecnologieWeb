% Esempio progetto regolatore per sistema carello
% 
% Controlli Automatici T
% Giuseppe Notarstefano
% Guido Carnevale
% 12/2023
%
% PARAMETRI FISICI CARRELLO
%
% massa carrello          = 0.5 kg
% costante elastica molla = 80 N/m
% coefficiente di attrito = 10 Ns/m
%
% SPECIFICHE
% 
% -) Errore a regime in risposta a un gradino w(t) = 1(t) e d(t) = 1(t) pari a 0.1
%
% -) Attenuazione di almeno 30dB per d(t)
%       con [omega_d_min, omega_d_MAX] = [0.001,0.1]
%
% -) Attenuazione di almeno 60dB per n(t)
%       con [omega_n_min, omega_n_MAX] = [10^3,10^4]
%
% -) S% <= 6%
% -) Ta,1 <= 0.5 s
%
% M_f > 45 gradi

clear all; close all; clc

%solo per visualizzione, pulsazione minima e massima
omega_plot_min = 1e-2;
omega_plot_max = 1e5;

%% Parametri

mm   = 0.5; % [kg] Massa
kk   = 80;  % [N/m^2] Coefficiente di attrito dinamico
bb   = 10;  % [N/m] Coefficiente di rigidezza della molla
x_1e = 0.5;

%% Coppia di equilibrio
x_2e = 0;
u_e  = kk*x_1e^2;
x_e  = [x_1e;x_2e];

%% Funzione di Trasferimento

A = [0 1;-2*kk*x_1e/mm -bb/mm];
B = [0;1/mm];
C = [1 0];
D = 0;

s  = tf('s');
GG = C*inv(s*eye(2) - A)*B + D;

%% Specifiche

% ampiezze gradini
WW = 1;
DD = 1;

% errore a regime
e_star = 0.1;

% attenuazione disturbo sull'uscita
A_d = 30;
omega_d_min = 0.001;
omega_d_MAX = 0.1;

% attenuazione disturbo di misura
A_n = 60;
omega_n_min = 1e4;
omega_n_MAX = 1e5;

% Sovraelongazione massima e tempo d'assestamento all'1%
S_star = 6;
T_star = 0.1;

% Margine di fase
Mf_star_esp = 45;


%% Diagramma di Bode

figure(1);
bode(GG,{omega_plot_min,omega_plot_max});
grid on, zoom on;

% return;


%% Regolatore statico - proporzionale senza poli nell'origine

% valore minimo prescritto per L(0)
mu_s_error = (DD+WW)/e_star;
mu_s_dist  = 10^(A_d/20);

% guadagno minimo del regolatore ottenuto come L(0)/G(0)
G_0 = abs(evalfr(GG,0));
G_omega_d_MAX = abs(evalfr(GG,j*omega_d_MAX));


RR_s = max(mu_s_error/G_0,mu_s_dist/G_omega_d_MAX); % RR_s = 800

% Sistema esteso
GG_e = RR_s*GG;


%% Diagrammi di Bode di Ge con specifiche
figure(2);
hold on;

% Mappatura sovraelongazione => Margine di fase
% Agire su xi_star per ridurre sovraelongazione
xi_star = abs(log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);
Mf_star = max(xi_star*100,Mf_star_esp);

% Mappatura tempo di assestamento => pulsazione critca
omega_c_min = 460/(Mf_star*T_star); % omega_c >= 460/(Mf*T^*) ~ 4.6
omega_c_max = omega_n_min;

% Specifiche su d
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

Bnd_Ta_x = [omega_plot_min; omega_c_min; omega_c_min; omega_plot_min];
Bnd_Ta_y = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;



phi_low = -270; % lower bound per il plot (teoricamente - infinito)

Bnd_Mf_x = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y = [Mf_star - 180; Mf_star - 180; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);


%% Design del regolatore dinamico


omega_c_star = 200;

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

R_high_frequency = 1/(1 + s/4e3);

RR_d = (1 + tau*s)/(1 + alpha * tau*s)*R_high_frequency;

RR = RR_s*RR_d;

LL = RR*GG; % funzione di anello

figure(3);
hold on;

% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);

% STOP qui per sistema con controllore dinamico + specifiche
if 0
    return;
end


%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);

% Risposta al gradino
figure(4);

T_simulation = 2*T_star;
[y_step,t_step] = step(WW*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = evalfr(WW*FF,0);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[LV*(1+S_star/100),LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento all'1%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.01),LV*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.01),LV*(1+0.01),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

%% Check disturbo in uscita

% Funzione di sensitività
SS = 1/(1+LL);
figure(5);

% Simulazione disturbo in uscita a pulsazione 0.05
omega_d = 0.05;
tt = 0:1e-2:2e2;
dd = DD*cos(omega_d*tt);
y_d = lsim(SS,dd,tt);
hold on, grid on, zoom on
plot(tt,dd,'m')
plot(tt,y_d,'b')
grid on
legend('d(t)','y_d(t)')

%% Check disturbo di misura

% Funzione di sensitività complementare
FF = LL/(1+LL);
figure(6);

% Simulazione disturbo di misura a pulsazione 1000
omega_n = omega_n_min;
NN      = 1;
tt = 0:1e-5:2*1e-3;
nn = NN*cos(omega_n*tt);
y_n = lsim(FF,nn,tt);
hold on, grid on, zoom on
plot(tt,nn,'m')
plot(tt,y_n,'b')
grid on
legend('n(t)','y_n(t')
