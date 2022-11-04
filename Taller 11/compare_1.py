import numpy as np
from propagator import propagator
import matplotlib.pyplot as plt
from drawnow import drawnow
# from propagator import Pt
import timeit
# define parameters

dt = 0.004 
tEnd = 3
Nt = 750
t = np.linspace(0, tEnd-dt, Nt)
Nz = 150 
Nx = 300
vp = 3000*np.ones((Nx,Nz))
dz = 20
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
stop=timeit.default_timer()
print('Time: ', stop - start)

# Load analytical solution and compare to numerical
Ptx = np.arange(20,280)
Pt =  np.transpose(np.reshape(P[Ptx, 100, :], [260, Nt]))

analitica = np.loadtxt( 'analytical.txt')
numerical = np.amax(analitica[:,1]) / np.amax (Pt[:,199]) * Pt[:,199]
fig, ax = plt.subplots()
ax.plot(analitica[:,0], analitica[:,1], label= 'Analítica', c= 'orange', ls= '--') 
ax.plot(analitica[:,0], numerical, label= 'Numérica', c= 'r')
ax.set_xlabel(f'Tiempo $[s]$')
ax.set_ylabel( f'Amplitud')
ax.legend()
plt.show()

# wave propagation movie

def make_fig():
    plt.imshow(P[:,:,i].T)
    plt.title('Iteration ' + str(i))

for i in range(0, P.shape[-1],5):
    drawnow(make_fig)
plt.show()