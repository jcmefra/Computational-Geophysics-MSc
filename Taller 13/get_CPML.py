
import numpy as np
from numba import jit

@jit
def get_CPML(CPMLimit,R,Vcpml,Nx,Nz,dx,dz,dt,frec):
    D_pml_x = CPMLimit*dx 
    D_pml_z = CPMLimit*dz
    d0_x = -3/(2*D_pml_x)*np.log(R)
    d0_z = -3/(2*D_pml_z)*np.log(R)
    x = np.zeros(CPMLimit+1)
    z = np.zeros(CPMLimit+1)
    alpha_x = np.zeros(CPMLimit+1)
    alpha_z = np.zeros(CPMLimit+1)
    x_half = np.zeros(CPMLimit+1)
    z_half = np.zeros(CPMLimit+1)
    alpha_x_half = np.zeros(CPMLimit+1)
    alpha_z_half = np.zeros(CPMLimit+1)
    for j in range(1,CPMLimit+2):
        x[j-1] = (CPMLimit-j+1)*dx
        z[j-1] = (CPMLimit-j+1)*dz
        alpha_x[j-1] = np.pi*frec*(D_pml_x-x[j-1])/D_pml_x
        alpha_z[j-1] = np.pi*frec*(D_pml_z-z[j-1])/D_pml_z
        x_half[j-1] = (CPMLimit-j+1)*dx-dx/2
        z_half[j-1] = (CPMLimit-j+1)*dz-dz/2
        alpha_x_half[j-1] = np.pi*frec*(D_pml_x-x_half[j-1])/D_pml_x
        alpha_z_half[j-1] = np.pi*frec*(D_pml_z-z_half[j-1])/D_pml_z
    #Bottom side
    d_z = np.zeros(Nz)
    b_z = np.zeros(Nz)
    a_z = np.zeros(Nz)
    d_z_half = np.zeros(Nz-1)
    b_z_half = np.zeros(Nz-1)
    a_z_half = np.zeros(Nz-1)
    for i in range(Nz-CPMLimit,Nz+1):
        d_z[i-1] = d0_z*Vcpml*((z[Nz-i]/D_pml_z)**2)
        b_z[i-1] = np.exp(-(d_z[i-1]+alpha_z[Nz-i])*dt)
        a_z[i-1] = d_z[i-1]/(d_z[i-1]+alpha_z[Nz-i])*(b_z[i-1]-1)
        if (i==Nz-CPMLimit):
            d_z_half[i-2] = 0
            b_z_half[i-2] = 0
            a_z_half[i-2] = 0
        else:
            d_z_half[i-2] = d0_z*(Vcpml)*((z_half[Nz-i]/D_pml_z)**2)
            b_z_half[i-2] = np.exp(-(d_z_half[i-2]+alpha_z_half[Nz-i])*dt)
            a_z_half[i-2] = d_z_half[i-2]/(d_z_half[i-2]+alpha_z_half[Nz-i])*(b_z_half[i-2]-1)
    #Left side
    d_x = np.zeros(Nx)
    b_x = np.zeros(Nx)
    a_x = np.zeros(Nx)
    d_x_half = np.zeros(Nx-1)
    b_x_half = np.zeros(Nx-1)
    a_x_half = np.zeros(Nx-1)
    for k in range(CPMLimit+1):
        d_x[k] = d0_x*Vcpml*((x[k]/D_pml_x)**2)
        b_x[k] = np.exp(-(d_x[k]+alpha_x[k])*dt)
        a_x[k] = d_x[k]/(d_x[k]+alpha_x[k])*(b_x[k]-1)
        if (k==CPMLimit):
            d_x_half[k] = 0
            b_x_half[k] = 0
            a_x_half[k] = 0
        else:
            d_x_half[k] = d0_x*(Vcpml)*((x_half[k]/D_pml_x)**2)
            b_x_half[k] = np.exp(-(d_x_half[k]+alpha_x_half[k])*dt)
            a_x_half[k] = d_x_half[k]/(d_x_half[k]+alpha_x_half[k])*(b_x_half[k]-1)
    #Right side
    for q in range(Nx-CPMLimit,Nx+1):
        d_x[q-1] = d0_x*Vcpml*((x[Nx-q]/D_pml_x)**2) 
        b_x[q-1] = np.exp(-(d_x[q-1]+alpha_x[Nx-q])*dt)
        a_x[q-1] = d_x[q-1]/(d_x[q-1]+alpha_x[Nx-q])*(b_x[q-1]-1)
        if (q==Nx-CPMLimit):
            d_x_half[q-2] = 0
            b_x_half[q-2] = 0
            a_x_half[q-2] = 0
        else:
            d_x_half[q-2] = d0_x*Vcpml*((x_half[Nx-q]/D_pml_x)**2) 
            b_x_half[q-2] = np.exp(-(d_x_half[q-2]+alpha_x_half[Nx-q])*dt)
            a_x_half[q-2] = d_x_half[q-2]/(d_x_half[q-2]+alpha_x_half[Nx-q])*(b_x_half[q-2]-1) 
    return a_x,a_x_half,b_x,b_x_half,a_z,a_z_half,b_z,b_z_half 
