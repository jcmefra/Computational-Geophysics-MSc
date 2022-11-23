function [Pt,P,d2P_dt2]=propagator(m, src, Ix0, Iz0, dx, dz, dt, max_offset, frec)

max_ix=max_offset/dx;
[Nx,Nz] = size( m );
Nt = size( src,1 );

% Initialize arrays
P = zeros( Nx,Nz,Nt ); P_Iz0 = zeros(Nt,Nx);
d2P_dx2 = 0; d2P_dz2 = 0;
d2P_dt2 = P ; dP_dx = d2P_dt2; dP_dz = d2P_dt2;
F_dPdx = d2P_dt2; F_dPdz = d2P_dt2; F_d2Pdx2 = d2P_dt2; F_d2Pdz2 = d2P_dt2;

v2 = m.^2;

l_att = 20; CPMLimit = l_att;
one_over_dx2 = 1/dx^2;
one_over_dz2 = 1/dz^2;
one_over_dx = 1/dx;
one_over_dz = 1/dz;

Vmax=max(max(m)); Vmin=min(min(m));
w_l_ = Vmin/frec; Points_per_wavelength = w_l_/dx
Courant_number = Vmax * dt/dx * sqrt(2)

% Functions for C-PML conditions
R = 1e-3;  Vcpml=Vmax;
[ a_x,a_x_half,b_x,b_x_half,a_z,a_z_half,b_z,b_z_half ] = get_CPML( CPMLimit,R,Vcpml,Nx,Nz,dx,dz,dt,frec );

for it = 2:Nt-1
    
    % Snapshot
    P_tmp=P(:,:,it);
    
    % Free-surface conditions
    %d2P_dt2( :,2,it ) = 0 ;
    P( :,2,it ) = 0 ;
    
    % Spatial finite differences
    for ix = CPMLimit+2:Nx-CPMLimit-1
        for iz = 3:Nz-CPMLimit-1
            d2P_dx2 = P_tmp( ix-1,iz ) - 2*P_tmp( ix,iz ) + P_tmp( ix+1,iz ) ; d2P_dx2 = one_over_dx2*d2P_dx2 ;
            d2P_dz2 = P_tmp( ix,iz-1 ) - 2*P_tmp( ix,iz ) + P_tmp( ix,iz+1 ) ; d2P_dz2 = one_over_dz2*d2P_dz2 ;
            
            d2P_dt2( ix,iz,it ) = d2P_dx2+d2P_dz2 ;
            
        end
    end
    P_tmp(Nx-1,:)=0; P_tmp(2,:)=0; P_tmp(:,Nz-1)=0;
    
    %% Begin: C-PML conditions
    for ix = 2:CPMLimit+1
        for iz = 2:Nz-1
            dP_dx( ix,iz ) = P_tmp( ix+1,iz ) - P_tmp( ix,iz ) ; dP_dx( ix,iz ) = one_over_dx*dP_dx( ix,iz ) ;
            dP_dz( ix,iz ) = P_tmp( ix,iz+1 ) - P_tmp( ix,iz ) ; dP_dz( ix,iz ) = one_over_dz*dP_dz( ix,iz ) ;
            
            % Here it should be half point not b_x
            F_dPdx(ix,iz) = F_dPdx(ix,iz)*b_x_half(ix) + a_x_half(ix)*dP_dx(ix,iz); dP_dx(ix,iz) = dP_dx(ix,iz)+F_dPdx(ix,iz);
            F_dPdz(ix,iz) = F_dPdz(ix,iz)*b_z(iz) + a_z(iz)*dP_dz(ix,iz); dP_dz(ix,iz) = dP_dz(ix,iz)+F_dPdz(ix,iz);
            
            d2P_dx2 = one_over_dx * ( dP_dx(ix,iz)-dP_dx(ix-1,iz) );
            d2P_dz2 = one_over_dz * ( dP_dz(ix,iz)-dP_dz(ix,iz-1) );
            
            F_d2Pdx2(ix,iz) = F_d2Pdx2(ix,iz)*b_x(ix) + a_x(ix)*d2P_dx2; d2P_dx2 = d2P_dx2+F_d2Pdx2(ix,iz);
            F_d2Pdz2(ix,iz) = F_d2Pdz2(ix,iz)*b_z(iz) + a_z(iz)*d2P_dz2; d2P_dz2 = d2P_dz2+F_d2Pdz2(ix,iz);
            
            d2P_dt2( ix,iz,it ) = d2P_dx2+d2P_dz2 ;
        end
    end
    for ix = Nx-1:-1:Nx-CPMLimit
        for iz = 2:Nz-1
            dP_dx( ix,iz ) = P_tmp( ix,iz ) - P_tmp( ix-1,iz ) ; dP_dx( ix,iz ) = one_over_dx*dP_dx( ix,iz ) ;
            dP_dz( ix,iz ) = P_tmp( ix,iz+1 ) - P_tmp( ix,iz ) ; dP_dz( ix,iz ) = one_over_dz*dP_dz( ix,iz ) ;
            
            % Here it should be half point not b_x
            F_dPdx(ix,iz) = F_dPdx(ix,iz)*b_x_half(ix-1) + a_x_half(ix-1)*dP_dx(ix,iz); dP_dx(ix,iz) = dP_dx(ix,iz)+F_dPdx(ix,iz);
            F_dPdz(ix,iz) = F_dPdz(ix,iz)*b_z(iz) + a_z(iz)*dP_dz(ix,iz); dP_dz(ix,iz) = dP_dz(ix,iz)+F_dPdz(ix,iz);
            
            d2P_dx2 = one_over_dx * ( dP_dx(ix+1,iz)-dP_dx(ix,iz) );
            d2P_dz2 = one_over_dz * ( dP_dz(ix,iz)-dP_dz(ix,iz-1) );
            
            F_d2Pdx2(ix,iz) = F_d2Pdx2(ix,iz)*b_x(ix) + a_x(ix)*d2P_dx2; d2P_dx2 = d2P_dx2+F_d2Pdx2(ix,iz);
            F_d2Pdz2(ix,iz) = F_d2Pdz2(ix,iz)*b_z(iz) + a_z(iz)*d2P_dz2; d2P_dz2 = d2P_dz2+F_d2Pdz2(ix,iz);
            
            d2P_dt2( ix,iz,it ) = d2P_dx2+d2P_dz2 ;
        end
    end
    for ix = CPMLimit+2:Nx-CPMLimit-1
        for iz = Nz-1:-1:Nz-CPMLimit
            dP_dx( ix,iz ) = P_tmp( ix+1,iz ) - P_tmp( ix,iz ) ; dP_dx( ix,iz ) = one_over_dx*dP_dx( ix,iz ) ;
            dP_dz( ix,iz ) = P_tmp( ix,iz ) - P_tmp( ix,iz-1 ) ; dP_dz( ix,iz ) = one_over_dz*dP_dz( ix,iz ) ;
            
            % Here it should be half point not b_x
            F_dPdx(ix,iz) = F_dPdx(ix,iz)*b_x(ix) + a_x(ix)*dP_dx(ix,iz); dP_dx(ix,iz) = dP_dx(ix,iz)+F_dPdx(ix,iz);
            F_dPdz(ix,iz) = F_dPdz(ix,iz)*b_z_half(iz-1) + a_z_half(iz-1)*dP_dz(ix,iz); dP_dz(ix,iz) = dP_dz(ix,iz)+F_dPdz(ix,iz);
            
            d2P_dx2 = one_over_dx * ( dP_dx(ix,iz)-dP_dx(ix-1,iz) );
            d2P_dz2 = one_over_dz * ( dP_dz(ix,iz+1)-dP_dz(ix,iz) );
            
            F_d2Pdx2(ix,iz) = F_d2Pdx2(ix,iz)*b_x(ix) + a_x(ix)*d2P_dx2; d2P_dx2 = d2P_dx2+F_d2Pdx2(ix,iz);
            F_d2Pdz2(ix,iz) = F_d2Pdz2(ix,iz)*b_z(iz) + a_z(iz)*d2P_dz2; d2P_dz2 = d2P_dz2+F_d2Pdz2(ix,iz);
            
            d2P_dt2( ix,iz,it ) = d2P_dx2+d2P_dz2 ;
        end
    end
    %% End: C-PML conditions
    
    d2P_dt2(:,:,it) = v2.*d2P_dt2(:,:,it);
    
    % Time integration
    P(:,:,it+1) = ( 2*P(:,:,it)-P(:,:,it-1) ) + dt*dt*d2P_dt2(:,:,it);
    
    % Source injection
    for is=Ix0
        P( is,Iz0,it+1 ) = P( is,Iz0,it+1 ) + src( it+1,is-min(Ix0)+1 );
    end
    
    % Surface component selection
    P_Iz0(it,:) = P(:,Iz0,it);
    
    % Counting solved percentage
    if mod(it/(Nt)*100,20) == 0
        co = num2str(it/(Nt)*100); display([co '%']);
    end
    
end

if (length(Ix0)>1)
    Pt = reshape(P(Ix0,Iz0,:),[length(Ix0),Nt])';
else
    if(Ix0>max_ix)
        Pt = reshape(P(max(21,Ix0-max_ix):min(Nx-20,Ix0+max_ix)-1,Iz0,:),[min(Nx-20,Ix0+max_ix)-max(21,Ix0-max_ix),Nt])';
    else
        Pt = reshape(P(max(21,Ix0-max_ix):min(Nx-20,Ix0+max_ix),Iz0,:),[min(Nx-20,Ix0+max_ix)-max(21,Ix0-max_ix)+1,Nt])';
    end
    
end
