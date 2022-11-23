function [ a_x,a_x_half,b_x,b_x_half,a_z,a_z_half,b_z,b_z_half ] = get_CPML( CPMLimit,R,Vcpml,Nx,Nz,dx,dz,dt,frec )

D_pml_x = CPMLimit*dx; D_pml_z = CPMLimit*dz;
d0_x = -3/(2*D_pml_x)*log(R);d0_z = -3/(2*D_pml_z)*log(R);
for j = 1:CPMLimit+1
    x(j) = (CPMLimit-j+1)*dx;
    z(j) = (CPMLimit-j+1)*dz;
    alpha_x(j) = pi*frec*(D_pml_x-x(j))/D_pml_x;
    alpha_z(j) = pi*frec*(D_pml_z-z(j))/D_pml_z;
        x_half(j) = (CPMLimit-j+1)*dx-dx/2;
        z_half(j) = (CPMLimit-j+1)*dz-dz/2;
        alpha_x_half(j) = pi*frec*(D_pml_x-x_half(j))/D_pml_x;
        alpha_z_half(j) = pi*frec*(D_pml_z-z_half(j))/D_pml_z;
end
% Bottom side
for j = Nz-CPMLimit:Nz
    d_z(j) = d0_z*Vcpml*((z(Nz-j+1)/D_pml_z)^2);
    b_z(j) = exp(-(d_z(j)+alpha_z(Nz-j+1))*dt);
    a_z(j) = d_z(j)/(d_z(j)+alpha_z(Nz-j+1))*(b_z(j)-1);
    if (j == Nz-CPMLimit)
       d_z_half(j-1) = 0;
       b_z_half(j-1) = 0;
       a_z_half(j-1) = 0;
    else
       d_z_half(j-1) = d0_z*(Vcpml)*((z_half(Nz-j+1)/D_pml_z)^2);
       b_z_half(j-1) = exp(-(d_z_half(j-1)+alpha_z_half(Nz-j+1))*dt);
       a_z_half(j-1) = d_z_half(j-1)/(d_z_half(j-1)+alpha_z_half(Nz-j+1))*(b_z_half(j-1)-1);
    end
end
% Left side
for i = 1:CPMLimit+1
    d_x(i) = d0_x*Vcpml*((x(i)/D_pml_x)^2);
    b_x(i) = exp(-(d_x(i)+alpha_x(i))*dt);
    a_x(i) = d_x(i)/(d_x(i)+alpha_x(i))*(b_x(i)-1);
    if (i == CPMLimit+1)
        d_x_half(i) = 0; 
        b_x_half(i) = 0;
        a_x_half(i) = 0;
    else
        d_x_half(i) = d0_x*(Vcpml)*((x_half(i)/D_pml_x)^2) ; 
        b_x_half(i) = exp(-(d_x_half(i)+alpha_x_half(i))*dt);
        a_x_half(i) = d_x_half(i)/(d_x_half(i)+alpha_x_half(i))*(b_x_half(i)-1);
    end
end
% Right side
for i = Nx-CPMLimit:Nx
    d_x(i) = d0_x*Vcpml*((x(Nx-i+1)/D_pml_x)^2);  
    b_x(i) = exp(-(d_x(i)+alpha_x(Nx-i+1))*dt);
    a_x(i) = d_x(i)/(d_x(i)+alpha_x(Nx-i+1))*(b_x(i)-1);
    if (i == Nx-CPMLimit)
       d_x_half(i-1) = 0;  
       b_x_half(i-1) = 0;
       a_x_half(i-1) = 0;
    else
       d_x_half(i-1) = d0_x*(Vcpml)*((x_half(Nx-i+1)/D_pml_x)^2)  ;
       b_x_half(i-1) = exp(-(d_x_half(i-1)+alpha_x_half(Nx-i+1))*dt);
       a_x_half(i-1) = d_x_half(i-1)/(d_x_half(i-1)+alpha_x_half(Nx-i+1))*(b_x_half(i-1)-1);
    end
end