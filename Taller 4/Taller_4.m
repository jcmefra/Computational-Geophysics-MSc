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
Rank_Gw = rank(Gw);
for i = 1:100
    P = 243-i;
    Up = U(:,1:P);
    Vp = V(:,1:P);
    Sp = S(1:P,1:P);
    G_dagger = Vp*inv(Sp)*Up';
    m_dagger (:,i) = G_dagger*dw; %G_dagger*dn;
    Nm(i) = norm(Vp*inv(Sp)*Up'*dw);
    Nresidual(i) = norm((Up*(Sp*Vp'*m_dagger(:,i)))-dw); 
end

%Matriz de resolución
[Ux,Sx,Vx] = svds(Gw,(243-15));
Gdagger_x = Vx*inv(Sx)*Ux';
Rm = Gdagger_x*Gw;
Traza = trace(Rm);

%Curva L
semilogx(Nresidual,Nm)
title("Curva L")
xlabel('Norma residual ||Gm - d||_2')
ylabel('Norma de la solución ||m||_2')
hold on
scatter(3.002401293390646e-06,0.005453781475499)
text (3.002401293390646e-06,0.005453781475499, '$\rightarrow P = 228$', 'Interpreter','latex')


%Slowness
figure, imagesc(reshape(m_dagger(:,15),16,16))
title('Matriz de slowness usando TSVD')
colormap bone
colorbar

figure, contourf(reshape(m_dagger(:,15),16,16))
title('inversión del bloque.')
colormap('bone')
colorbar

min = min(m_dagger(:,15));
max = max(m_dagger(:,15));

figure
mesh (Rm)
title('Matriz de resolución')
colorbar

