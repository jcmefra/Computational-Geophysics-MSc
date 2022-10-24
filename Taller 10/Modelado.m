%Computational Geophysics Course
%Dr.Ing. Sergio Abreo.
%Acustical Wave equation implementation using constant density.


%Finite diferences in time domain discretization.
%Pseudocode discussed in class
%Data input.

Tout=1000; %Propagation time
Nx=200;   %Number of point in the x direction
Nz=200;   %Number of point in the z direction
c=1200.0; %Background velocity.
dh=2.0;   %spatial step
dt=dh/(c*sqrt(2));
G=c*dt/dh; %Courant parameter.

P1 = zeros(Nz,Nx);
P2 = zeros(Nz,Nx);
P3 = zeros(Nz,Nx);

alpha=0.15;       %alpha defines the pulse width

x0 = dh*Nx/2;
z0 = dh*Nz/3;

video = zeros(Nz,Nx,Tout);

for i=1:Nz
    for j=1:Nx
        r2 = (dh*i-z0)*(dh*i-z0)+(dh*j-x0)*(dh*j-x0);
        P2(i,j) = sin(1-alpha*r2)*exp(-alpha*r2);    % Source
    end
end

for t=1:Tout
    for i=2:Nz-1
        for j=2:Nx-1
            P3(i,j) = (2-4*G*G)*P2(i,j)+G*G*(P2(i+1,j)+P2(i-1,j)+P2(i,j+1)+P2(i,j-1))-P1(i,j);
        end
    end
    
    for i=2:Nz-1
        for j=2:Nx-1
            P1(i,j) = P2(i,j);
            P2(i,j) = P3(i,j);
        end
    end
    
    for i=1:Nz
        for j=1:Nx
            video(i,j,t) = P3(i,j);
        end
    end
end

paso=1;       % Video frames step.
figure

for i=1:paso:Tout
    imagesc(video(:,:,i)),colorbar % Reading the forward modeling
    hold on
    xlabel('Distance (m)');
    ylabel('Depth (m)');
    str=sprintf('Iteration %d',i);title(str);
    pause(0.001)
end