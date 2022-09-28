close all
clear all
clc

load("jitorres_crosswell.mat")
sigma = 0.5;
for i = 1:256
    v(i) = 1/sigma;
end
%Matriz
W = diag(v,0);
Gw = W*G;
dw = W*dn;
[U,S,V] = svd(Gw);
Alpha = 0;
%Ciclo Tikhonov
for j = 1:50
    Alpha = (j*0.1)-0.009;
    m_a = 0;
    for i = 1:256
        si = S(i,i);
        fi = (si.^2)./((si.^2)+Alpha^2);
        Ui = U(:,i);
        Vi = V(:,i);
        m_a = m_a+(fi*Ui'*dw*Vi)/si;
    end
    m_alpha(:,j) = m_a;
    Nm_alpha(j)=norm(m_alpha(:,j));
    Nresidual_alpha(j)=norm(Gw*m_alpha(:,j)-dw);

end

%Matriz de resolución 
G_sharpx = inv(Gw' * Gw + 0.5^2 * eye(256))*Gw';
Rm = G_sharpx*Gw;
Traza = trace(Rm);

%Plot L-curve
semilogx(Nresidual_alpha,Nm_alpha)
title('Curva L (con alfa)')
xlabel('Norma residual ||Gm-d||_{2}')
ylabel('Norma de la solución ||m||_2')
hold on
scatter(1.14435e-05,0.00545335)
text (1.14435e-05,0.00545335, '$\rightarrow \alpha = 0.5$', 'Interpreter','latex')

%Plot Slowness
figure, imagesc(reshape(m_alpha(:,5),16,16))
title('inversión del bloque. Min: 2.9723e-4. Max; 3.5680e-4')
colormap bone
colorbar

figure, contourf(reshape(m_alpha(:,5),16,16))
title('inversión del bloque.')
colormap bone
colorbar

min = min(m_alpha(:,5));
max = max(m_alpha(:,5));

%Plot res matrix
figure
mesh (Rm)
title('Matriz de resolución')
colorbar
