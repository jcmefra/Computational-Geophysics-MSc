close all
clear all
clc
load ('jitorres_crosswell.mat')
sigma=0.5;
for i=1:256
v(i)=1/sigma ;
end
%Matriz
W=diag(v,0);
Gw=W*G;
dw=W*dn;
%Algoritmo SIRT
%G1=Gw>0;
G1=Gw~=0;
m=zeros(256,1);
N=zeros(256,1);
L=zeros(256,1);
tic
for i=1:256
    N(i)=sum(G1(i,:));
    L(i)=sum(Gw(i,:));
end
res=1;
l=100;
K=sum(G1,1);
Delta=zeros(256,1);
for z=1:500
for i=1:256
    q=G1(i,:)*m;
    tiempos=dw(i)/L(i) - q/(l*N(i));
    Delta=Delta+tiempos*G1(i,:)';
end
m=m+Delta./K';
res(z)=norm(G1*m-dw);
res1(z)=norm(m,1);
end
m_escalado = m/100;
toc
figure, imagesc(reshape(m_escalado,16,16))
title('Matriz de slowness usando SIRT')
colormap bone
colorbar
figure, contourf(reshape(m_escalado,16,16))
title('Mapa de contornos de slowness usando SIRT')
colormap bone
colorbar
figure, semilogx(res)
title('Norma de los residuales por iteración')
xlabel('Iteraciones')
ylabel('Norma Residual ||Gm-d||_2')
figure, loglog(res1)
title('Norma de la solución por iteración')
xlabel('Iteraciones')
ylabel('Norma ||m||_2')