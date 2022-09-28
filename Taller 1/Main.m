x = [6; 10.1333; 14.2667; 18.4; 22.5333; 26.6667];
t = [3.4935; 4.2853; 5.1374; 5.8181; 6.8632; 8.1841];

G = [ones(length(x),1) x];

% Solución paso por paso
A = G'*G;
B = inv(A);  % Matriz de covarianza
% Cálculo de parámetros
m = B*G'*t;

% Calculo 95th percentil de chi cuadrado
DELTA2 = chi2inv(0.95, 2);
disp(DELTA2)

plot_ellipse(DELTA2, B, m, 10000)