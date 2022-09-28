% Método de míminos cuadrados Ax = b. Donde A es G y b es d
G = [1 6; 1 10.1333; 1 14.2667; 1 18.4; 1 22.5333; 1 26.6667];
x = [6; 10.1333; 14.2667; 18.4; 22.5333; 26.6667];
t = [3.4935; 4.2853; 5.1374; 5.8181; 6.8632; 8.1841];

[m,flag,relres,iter,resvec,lsvec] = lsqr(G,t);

% Gráfico del residual de mínimos cuadrados
N = length(resvec);
% semilogy(0:N-2,lsvec,'--o',0:N-1,resvec,'-o')
% legend("Least-squares residual","Relative residual")

% Solución paso por paso
Gt = transpose(G);

ml = inv(Gt * G)*Gt*t;

% Grafico de los datos
subplot(1,2,1)
scatter(x,t)
hold on
xlabel('Distancia [km]')
ylabel('Tiempo [seg]')
title('Regresión lineal experimento de refracción sísmica')
grid on

% Grafico de la línea de mejor ajuste
yCalc = G*ml;
plot(x,yCalc,'--')
legend('Datos','Mejor ajuste','Location','southeast');

% Residuales
subplot(1,2,2)
r = t - yCalc;
scatter(t, r, 'red', 'filled')
xlabel('Tiempo medido [s]')
ylabel('Valor residual (d-Gm)')
legend('Residuales', 'Location', 'southeast')
title('Valores residuales vs Tiempo medido')
grid on

