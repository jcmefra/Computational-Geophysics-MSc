close all
clear all
clc

%Laplaciano
L = zeros (14*14,256);
k = 1;
for i = 2:15
    for j = 2:15
        N = zeros (16,16);
        N(i,j) = -4;
        N(i,j+1) = 1;
        N(i,j-1) = 1;
        N(i+1,j) = 1;
        N(i-1,j) = 1;
        L(k,:) = (reshape(N,256,1))';
        k = k+1;
    end
end

load ("jitorres_crosswell.mat")
sigma = 0.5;
for i = 1:256
    v(i) = 1/sigma;
end
W = diag(v,0);
Gw = W*G;
dw = W*dn;
[U,V,X,C,M] = gsvd(Gw,L);
Alpha = 0.1;

m_soluciones = [];
elementos_plot = [];

for j = 0.1:0.1:10
    alpha = j;

    G_sharp = inv(Gw' * Gw + (alpha^2)*((L')*L))*Gw';
    
    mL_alpha = G_sharp*dw;
    Norma_soluciones = norm(L*mL_alpha);
    NmL_residual = Gw*mL_alpha - dw;
    Norma_residual = norm(NmL_residual);
    m_soluciones = [m_soluciones mL_alpha];

    a = [alpha Norma_residual Norma_soluciones];
    elementos_plot = [elementos_plot; a];
end

G_sharpx = inv(Gw' * Gw + (0.5^2 * (L')*L))*Gw';
Rm = G_sharpx*Gw;
Traza = trace(Rm);
M = max(Rm, [], 'all');

figure
semilogx(elementos_plot(:,2), elementos_plot(:,3))
title('Curva L (con alfa y L)')
xlabel('Norma residual ||Gm-d||_{2}')
ylabel('Norma de la solución ||Lm||_{2}')
hold on
scatter(elementos_plot(5,2), elementos_plot(5,3))
text (elementos_plot(5,2), elementos_plot(5,3), '$\rightarrow \alpha = 0.5$', 'Interpreter','latex')

figure, imagesc(reshape(m_soluciones(:,5),16,16))
title('inversión del bloque.')
colormap bone
colorbar

figure, contourf(reshape(m_soluciones(:,5),16,16))
title('inversión del bloque.')
colormap bone
colorbar

figure
mesh (Rm)
title('Matriz de resolución')
colorbar