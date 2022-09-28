% Clase minimos
x=[6; 10.1333; 14.2667; 18.4; 22.5333; 26.6667 ];
t=[3.4935; 4.2853; 5.1374; 5.8181; 6.8632; 8.1841];
% Para crear una matrix
G=[ones(6,1) x];
A=G'*G % La transpuesta de una matriz es G' GTGM=GTd
B=inv(A) %Matriz inversa de A usar inv ----  M=(GT*G)^1*GTd
M=B*G'*t  % Solucion de minimos cuadrados L2, una matriz que contiene la pendiente y el intercepto

scatter(x,t) % Graficar datos discretos
grid on
hold on

ti= M(1)+M(2)*x; %Graficar los puntos de la funcion ti
scatter(x,ti) % scatter es para discretos y el problema es discreto

t3=t-ti; % residuales
scatter(x,t3)
ylabel('Tiempos')
xlabel('Distancia')
title('Minimos Cuadrados')
legend('obs','pred','res')

%t3. es para tomar todos los puntos de la matriz
c=sum(t3.^2)

%p-value, chi2cdf(x, grados de libertad)
p=1-chi2cdf(c,4)

x = 0:0.2:18;
y = chi2pdf(x,4);

figure;
plot(x,y)
hold on
xl.LabelVerticalAlignment = 'bottom';
xline(c, '--red', 'Chi Cuadrado Obs')
xlabel('Observation')
ylabel('Probability Density')

