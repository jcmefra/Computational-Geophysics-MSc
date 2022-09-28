m = 22;
b=12;
x=linspace(0,25,30);
y=m*x + b;

%crear números aleatorios
n=30;
R = [0 30];
z = rand(n,1)*range(R)+min(R);
z = z';
y_error = y+z;

%construcción de outliers

m=100;
y_error(5) = y_error(5)+m;
y_error(16) = y_error(16)+m;
y_error(11) = y_error(11)+m;
y_error(23) = y_error(29)+m;
y_error(29) = y_error(29)+m;

%Matriz de unos, valores de X
G=[ones(30,1) x'];

A=G'*G; % La transpuesta de una matriz es G' GTGM=GTd
B=inv(A); %Matriz inversa de A usar inv ----  M=(GT*G)^1*GTd
M=B*G'*(y_error)';  % Solucion de minimos cuadrados L2, una matriz que contiene la pendiente y el intercepto

yl2 = M(1)+M(2)*x'; %gráfica norma l2
r = y_error' - yl2; %residuales

%Norma L1_1 
R_1 = diag(1./(abs(r)))

ml1_1 = inv(G'*R_1*G)*G'*R_1*y_error';
yl1_1 = G*ml1_1;

r_2 =  y_error' - yl1_1;

%norma L1_2
R_2 = diag(1./(abs(r_2)))

ml1_2 = inv(G'*R_2*G)*G'*R_2*y_error';
yl1_2 = G*ml1_2;

r_3 =  y_error' - yl1_2;

%norma L1_3

R_3 = diag(1./(abs(r_3)))

ml1_3 = inv(G'*R_3*G)*G'*R_3*y_error';
yl1_3 = G*ml1_3;


scatter(x,y_error,'filled')
hold on
grid on
plot(x,y, '--k')
plot(x,yl2,'color', [0.8500 0.3250 0.0980])
plot(x,yl1_1, 'color', [0.6350 0.0780 0.1840])
plot(x,yl1_2, 'b')
plot(x,yl1_3, 'color', [0.4940 0.1840 0.5560])
ylabel('Y')
xlabel('X')
title('Normas L1 y L2')
legend('Valores','Recta ideal','Norma L2', 'Norma L1_1', 'Norma L1_2', 'Norma L1_3')