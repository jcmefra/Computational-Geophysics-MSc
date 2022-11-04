import numpy as np
from numba import jit
from get_CPML import get_CPML

@jit
def propagator(m, src, Ix0, Iz0, dx, dz, dt, max_offset, frec):
    max_ix = max_offset/dx
    Nx, Nz = m.shape
    Nt = len(src)
    
    # initialize arrays
    P = np.zeros((Nx, Nz, Nt))
    P_tmp = np.zeros((Nx, Nz))
    P_Iz0 = np.zeros((Nt, Nx))
    d2P_dx2 = 0
    d2P_dz2 = 0
    d2P_dt2 = np.zeros((Nx, Nz, Nt))
    dP_dx = np.zeros((Nx, Nz))
    dP_dz = np.zeros((Nx, Nz))
    F_dPdx = np.zeros((Nx, Nz))
    F_dPdz = np.zeros((Nx, Nz))
    F_d2Pdx2 = np.zeros((Nx, Nz))
    F_d2Pdz2 = np.zeros((Nx, Nz))
    v2 = m**2
    
    l_att = 20
    CPMLimit = l_att
    one_over_dx2 = 1/dx**2
    one_over_dz2 = 1/dz**2
    one_over_dx = 1/dx
    one_over_dz = 1/dz
    Vmax = np.max(m)
    Vmin = np.min(m)
    w_l_ = Vmin/frec 
    Points_per_wavelength = w_l_/dx
    Courant_number = Vmax*dt/dx*np.sqrt(2)
    #print(Courant_number)
    
    # Functions for C-PML conditions
    R = 1e-3
    Vcpml = Vmax
    a_x, a_x_half, b_x, b_x_half, a_z, a_z_half, b_z, b_z_half = get_CPML(CPMLimit,R,Vcpml,Nx,Nz,dx,dz,dt,frec)
    
    for it in range(1, Nt-1):
        # Snapshot
        P_tmp[:,:] = P[:,:,it]
        
        # Free-surface conditions
        d2P_dt2[:,1,it] = 0
        
        # spatial finite differences
        for ix in range(CPMLimit+2-1, Nx-CPMLimit-1):
            for iz in range(2, Nz-CPMLimit-1):
                d2P_dx2 = (P_tmp[ix+1,iz] + P_tmp[ix-1,iz] - 2*P_tmp[ix,iz])*one_over_dx2
                d2P_dz2 = (P_tmp[ix,iz+1] + P_tmp[ix,iz-1] - 2*P_tmp[ix,iz])*one_over_dz2
                d2P_dt2[ix,iz,it] = d2P_dx2 + d2P_dz2
        
        P_tmp[Nx-2,:] = 0
        P_tmp[1,:] = 0
        P_tmp[:,Nz-2] = 0
        
        # Begin: C-PML conditions
        for ix in range(1, CPMLimit+1):
            for iz in range(1, Nz-1):
                dP_dx[ix,iz] = P_tmp[ix+1,iz] - P_tmp[ix,iz]
                dP_dx[ix,iz] = one_over_dx*dP_dx[ix,iz]
                dP_dz[ix,iz] = P_tmp[ix,iz+1] - P_tmp[ix,iz]
                dP_dz[ix,iz] = one_over_dz*dP_dz[ix,iz]
            
                # Here it should be half point not b_x
                F_dPdx[ix,iz] = F_dPdx[ix,iz]*b_x_half[ix] + a_x_half[ix]*dP_dx[ix,iz]
                dP_dx[ix,iz] = dP_dx[ix,iz] + F_dPdx[ix,iz]
                F_dPdz[ix,iz] = F_dPdz[ix,iz]*b_z[iz] + a_z[iz]*dP_dz[ix,iz] 
                dP_dz[ix,iz] = dP_dz[ix,iz] + F_dPdz[ix,iz]
            
                d2P_dx2 = one_over_dx*(dP_dx[ix,iz] - dP_dx[ix-1,iz])
                d2P_dz2 = one_over_dz*(dP_dz[ix,iz] - dP_dz[ix,iz-1])
            
                F_d2Pdx2[ix,iz] = F_d2Pdx2[ix,iz]*b_x[ix] + a_x[ix]*d2P_dx2
                d2P_dx2 = d2P_dx2 + F_d2Pdx2[ix,iz]
                F_d2Pdz2[ix,iz] = F_d2Pdz2[ix,iz]*b_z[iz] + a_z[iz]*d2P_dz2
                d2P_dz2 = d2P_dz2 + F_d2Pdz2[ix,iz]
                d2P_dt2[ix,iz,it] = d2P_dx2 + d2P_dz2
            
        for ix in range(Nx-2, Nx-CPMLimit-2, -1):
            for iz in range(1, Nz-1):
                dP_dx[ix,iz] = P_tmp[ix,iz] - P_tmp[ix-1,iz]
                dP_dx[ix,iz] = one_over_dx*dP_dx[ix,iz]
                dP_dz[ix,iz] = P_tmp[ix,iz+1] - P_tmp[ix,iz]
                dP_dz[ix,iz] = one_over_dz*dP_dz[ix,iz]
            
                # Here it should be half point not b_x
                F_dPdx[ix,iz] = F_dPdx[ix,iz]*b_x_half[ix-1] + a_x_half[ix-1]*dP_dx[ix,iz]
                dP_dx[ix,iz] = dP_dx[ix,iz] + F_dPdx[ix,iz]
                F_dPdz[ix,iz] = F_dPdz[ix,iz]*b_z[iz] + a_z[iz]*dP_dz[ix,iz]
                dP_dz[ix,iz] = dP_dz[ix,iz] + F_dPdz[ix,iz]
            
                d2P_dx2 = one_over_dx*(dP_dx[ix+1,iz] - dP_dx[ix,iz])
                d2P_dz2 = one_over_dz*(dP_dz[ix,iz] - dP_dz[ix,iz-1])
          
                F_d2Pdx2[ix,iz] = F_d2Pdx2[ix,iz]*b_x[ix] + a_x[ix]*d2P_dx2
                d2P_dx2 = d2P_dx2 + F_d2Pdx2[ix,iz]
                F_d2Pdz2[ix,iz] = F_d2Pdz2[ix,iz]*b_z[iz] + a_z[iz]*d2P_dz2
                d2P_dz2 = d2P_dz2 + F_d2Pdz2[ix,iz]
            
                d2P_dt2[ix,iz,it] = d2P_dx2 + d2P_dz2
           
        for ix in range(CPMLimit+1, Nx-CPMLimit-1):
            for iz in range(Nz-2, Nz-CPMLimit-2, -1):
                dP_dx[ix,iz] = P_tmp[ix+1,iz] - P_tmp[ix,iz]
                dP_dx[ix,iz] = one_over_dx*dP_dx[ix,iz]
                dP_dz[ix,iz] = P_tmp[ix,iz] - P_tmp[ix,iz-1]
                dP_dz[ix,iz] = one_over_dz*dP_dz[ix,iz]
            
                # Here it should be half point not b_x
                F_dPdx[ix,iz] = F_dPdx[ix,iz]*b_x[ix] + a_x[ix]*dP_dx[ix,iz]
                dP_dx[ix,iz] = dP_dx[ix,iz] + F_dPdx[ix,iz]
                F_dPdz[ix,iz] = F_dPdz[ix,iz]*b_z_half[iz-1] + a_z_half[iz-1]*dP_dz[ix,iz]
                dP_dz[ix,iz] = dP_dz[ix,iz] + F_dPdz[ix,iz]
            
                d2P_dx2 = one_over_dx*(dP_dx[ix,iz] - dP_dx[ix-1,iz])
                d2P_dz2 = one_over_dz*(dP_dz[ix,iz+1] - dP_dz[ix,iz])
            
                F_d2Pdx2[ix,iz] = F_d2Pdx2[ix,iz]*b_x[ix] + a_x[ix]*d2P_dx2
                d2P_dx2 = d2P_dx2 + F_d2Pdx2[ix,iz]
                F_d2Pdz2[ix,iz] = F_d2Pdz2[ix,iz]*b_z[iz] + a_z[iz]*d2P_dz2
                d2P_dz2 = d2P_dz2 + F_d2Pdz2[ix,iz]
            
                d2P_dt2[ix,iz,it] = d2P_dx2 + d2P_dz2
            
        # End: C-PML conditions
    
        # Time integration
        P[:,:,it+1] = dt**2*v2*d2P_dt2[:,:,it] + 2*P[:,:,it] - P[:,:,it-1]
        # Source injection
        P[Ix0,Iz0,it+1] = P[Ix0,Iz0,it+1] + src[it+1]
        #print(src)
        # Surface component selection
        P_Iz0[it,:] = P[:,Iz0,it]
    
        # Counting solved percentage
        #if ((it/Nt)*100)%20 == 0:
           # co = str(it/(Nt)*100)
           # print(co + '%')

    o = np.arange(np.max(np.array([21, Ix0-max_ix])), np.min(np.array([Nx-20, Ix0+max_ix]))+1)
    o = o.astype(int)
    o2 = int(np.min(np.array([Nx-20, Ix0+max_ix])) - np.max(np.array([21, Ix0-max_ix])) + 1)

    Pt = np.reshape(P[o, Iz0, :], (o2, Nt), order="F").T
    
    return Pt, P, d2P_dt2