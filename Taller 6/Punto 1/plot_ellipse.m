function plot_ellipse(DELTA2,C,m,n)
% This function receives DELTA2 as the 95th percentil of a chi-square
% distribution
% C as the covariance matriz of your parameters
% m as the parameters found by the minimization of 2-norm
% n as the number of points on the ellipse

    %construct a vector of n equally-spaced angles from (0,2*pi)
    theta=linspace(0,2*pi,n)';

    %corresponding unit vector
    xhat=[cos(theta),sin(theta)];

    Cinv=inv(C);%inversa covarianza, Diagonalizando se encuentra la dirección de los semiejes
    %preallocate output array

    r=zeros(n,2);
    for i=1:n
        %store each (x,y) pair on the confidence ellipse
        %in the corresponding row of r
        r(i,:)=sqrt(DELTA2/(xhat(i,:)*Cinv*xhat(i,:)'))*xhat(i,:);
    end

    plot(m(1)+r(:,1), m(2)+r(:,2), 'black', 'LineWidth', 1.5);
    hold on
    plot(m(1),m(2),'k.','MarkerSize',14)
    hold on
    text(m(1),m(2), ' m^{\ast}','FontSize',15)
    hold on
    title('Elipsoide de error')
    xlabel('Coeficiente de almacenamiento S')
    ylabel('Transmisividad T')
    grid on
%   axis equal
end