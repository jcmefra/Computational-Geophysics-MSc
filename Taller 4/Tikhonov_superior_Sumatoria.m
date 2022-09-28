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
Alpha = 0;
for j = 1:50
    Alpha = (j*0.1)-0.009;
    m_a = 0;
    for i = 1:196
        mu = sqrt(diag(M'*M));
        lambda = sqrt(diag(C'*C));
        gamma = lambda./mu;
        fi = (gamma(i).^2)./((gamma(i).^2)+Alpha^2);
        Ui = U(:,i);
        Y = inv(X');
        Yi = Y(:,i);
        m_a = m_a + (fi*Ui'*dw*Yi)./lambda(i);
    end
    m_alpha(:,j) = m_a;
    NLm_alpha(j)=norm(L*m_alpha(:,j));
    Nresidual_alpha(j)=norm((Gw*m_alpha(:,j))-dw);
end

G_sharpx = inv(Gw' * Gw + (0.5^2 * (L')*L))*Gw';
Rm = G_sharpx*Gw;
Traza = trace(Rm);

%Plot L-curve
loglog(Nresidual_alpha,NLm_alpha)
title('Curva L (con alfa y L)')
xlabel('Norma residual ||Gm-d||_{2}')
ylabel('Norma de la solución ||Lm||_{2}')

%Plot Slowness
figure, imagesc(reshape(m_alpha(:,5),16,16))
title('inversión del bloque.')
colormap('default')
colorbar

min = min(m_alpha(:,5));
max = max(m_alpha(:,5));

%Plot res matrix
figure
mesh (Rm)
title('Matriz de resolución')
colorbar
