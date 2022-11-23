import numpy as np
from propagator import propagator
import matplotlib.pyplot as plt
from drawnow import drawnow
import timeit

# define parameters

dt = 0.004 
tEnd = 3
Nt = 750
t = np.linspace(0, tEnd-dt, Nt)
Nz = 150 
Nx = 300
vp = 3000*np.ones((Nx,Nz))
dz = 25
dx = dz
Ix0 = 100
Iz0 = 3
frec = 5

# source

a = (np.pi*5)**2
t0 = 1
g = (-2*a*(t-t0)*np.exp(-a*(t-t0)**2)).T

# compute numeral wave propagation

start=timeit.default_timer()
Pt, P, d2P_dt2 = propagator(vp, g, Ix0, Iz0, dx, dz, dt, 10e3, frec)

# # Load analytical solution and compare to numerical
# Ptx = np.arange(20,280)
# Pt =  np.transpose(np.reshape(P[Ptx, 100, :], [260, Nt]))

# analitica = np.loadtxt( 'analytical.txt')
# numerical = np.amax(analitica[:,1]) / np.amax (Pt[:,199]) * Pt[:,199]
# fig, ax = plt.subplots()
# ax.plot(analitica[:,0], analitica[:,1], label= 'Analítica', c= 'orange', ls= '--') 
# ax.plot(analitica[:,0], numerical, label= 'Numérica', c= 'r')
# ax.set_xlabel(f'Tiempo $[s]$')
# ax.set_ylabel( f'Amplitud')
# ax.legend()
# plt.show()

# shotgather
Sx1=25;
Sx2=100;
Sx3=150;
Sx4=200;
Sx5=275;

Pmod1, P, d2P_dt2 = propagator(vp, g, Sx1, Iz0, dx, dz, dt, 10e3, frec)
Pmod2, P, d2P_dt2 = propagator(vp, g, Sx2, Iz0, dx, dz, dt, 10e3, frec)
Pmod3, P, d2P_dt2 = propagator(vp, g, Sx3, Iz0, dx, dz, dt, 10e3, frec)
Pmod4, P, d2P_dt2 = propagator(vp, g, Sx4, Iz0, dx, dz, dt, 10e3, frec)
Pmod5, P, d2P_dt2 = propagator(vp, g, Sx5, Iz0, dx, dz, dt, 10e3, frec)

stop=timeit.default_timer()
print('Time: ', stop - start)

# fig, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1,5)
# fig.suptitle('Shot Gathers')
# ax1.imshow(Pmod1[:,:], vmin=-0.1, vmax=0.1)
# ax1.set_title('Shot gather 1')
# ax2.imshow(Pmod2[:,:], vmin=-0.1, vmax=0.1)
# ax2.set_title('Shot gather 2')
# ax3.imshow(Pmod3[:,:], vmin=-0.1, vmax=0.1)
# ax3.set_title('Shot gather 3')
# ax4.imshow(Pmod4[:,:], vmin=-0.1, vmax=0.1)
# ax4.set_title('Shot gather 4')
# ax5.imshow(Pmod5[:,:], vmin=-0.1, vmax=0.1)
# ax5.set_title('Shot gather 5')
# fig.subplots_adjust(wspace=1)
# plt.show()

# # wave propagation and dPdT movie
# fig, (ax1,ax2) = plt.subplots(1,2)
# for i in range(0, P.shape[-1], 25):
#     fig.suptitle('Iteration'+ str(i))
#     uno = ax1.imshow(P[:,:,i].T,vmin=-1, vmax=1)
#     ax1.set_title('P Wave')
#     dos = ax2.imshow(d2P_dt2[:,:,i].T,vmin=-0.001, vmax=0.001)
#     ax2.set_title('d2P_dt2 Wave')
#     plt.pause(0.001)
# plt.show()


# wave propagation movie
def make_fig():
    plt.imshow(P[:,:,i].T)
    plt.title('Iteration ' + str(i))

for i in range(0, P.shape[-1], 5):
    drawnow(make_fig)
plt.show()