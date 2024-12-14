%% FUNZIONE DI TRASFERIMENTO G con ZONE PROIBITE
s=tf('s');
G=0.001155/(s^2+0.00075*s+0.06454);
figure
hold on;
omega_plot_min = 1e-1;
omega_plot_max = 1e3;
degree_plot_min = -360;
degree_plot_max = 360;
dB_plot_min = -200;
db_plot_max = 200;
%omega_min = 1;
T_a_5_str = 0.03;
S_str = 0.02;
mu_min = 399;
xi_min = (abs(log(S_str)))/(sqrt(pi^2+log(S_str)^2))+0.2;
omega_c_min = 3/(xi_min*T_a_5_str);
M_f_min = 30;

omega_min_d = omega_plot_min;
omega_max_d = 0.5;
omega_min_n = 1e5;
omega_max_n = 5e6;
A_d = 40;
A_n = 63;


x_wc = [omega_plot_min; omega_plot_min; omega_c_min; omega_c_min];
y_wc = [dB_plot_min; 0; 0; dB_plot_min];
x_d = [omega_min_d; omega_min_d; omega_max_d; omega_max_d];
y_d = [A_d; 0; 0; A_d];
x_n = [omega_min_n; omega_min_n; omega_max_n; omega_max_n];
y_n = [-A_n; db_plot_max; db_plot_max; -A_n];
x_Mf = [omega_plot_min; omega_plot_min; omega_c_min; omega_c_min ];
y_Mf = [-M_f_min; degree_plot_min; degree_plot_min; -M_f_min];
patch(x_wc, y_wc, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_d, y_d, 'r', 'FaceAlpha', 0.1);
hold on;
patch(x_n, y_n, 'r', 'FaceAlpha', 0.1);
hold on;

bode(G);

hold on;
patch(x_Mf,y_Mf, 'r', 'FaceAlpha', 0.1);