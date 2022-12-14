from tqdm import tqdm
#import heuristics
from heuristics import punto1 
import numpy as np
import matplotlib.pyplot as plt

S_0 = []
T_0 = []
muestreo_plot = []
w = 0
for w in tqdm(range(0,200)): #200 essays m_0
    w+=1
    #xx=[0.002, 0.585]
    #xx=np.array(xx)
    #punto1(xx)
    intervalo1=0.05  # 32 ackl, 100 boha1,   5 drop, 10 matyas, 1 punto1, 10 punto2
    paso1=2       # 1 ackl,   1 boha1, 0.5 drop,  1 matyas, 0.1 punto1, 0.5 punto2
    intervalo2=1
    paso2=2

    # X,Y = np.meshgrid(np.arange(-intervalo1,intervalo1,paso1), np.arange(-intervalo2,intervalo2,paso2))
    # k1,k2= X.shape
    # Z=np.zeros((k1,k2))

    # for i in range(0, k1):
    #     for j in range(0, k2):
    #         a= np.array([X[i, j], Y[i,j]])
    #         Z[i][j]= punto1(a)
            
    T=10 
    temp1=(np.random.rand())*intervalo1*paso1
    temp2=(np.random.rand())*intervalo2*paso2
    m0=np.array([temp1, temp2])
    m0=m0.ravel()
    Em0=punto1(m0)

    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
    #view(-17,18)#ackley
    #view(-69,67)#boha1
    #view(-161,29)#drop
    #ax.view_init(30,angle)#matya

    # plt.plot(m0[0],m0[1],Em0,'ro',label='Punto inicial')

    while T > 0: # loop over the temperature T
        for i in range (100): # loop over a number of random moves/temperature
            temp1=(np.random.rand())*intervalo1*paso1
            temp2=(np.random.rand())*intervalo2*paso1
            m1=np.array([temp1, temp2])
            m1=m1.ravel()
            Em1=punto1(m1)
            Delta_e=Em1-Em0
            P=np.exp(-Delta_e/T)
            if Delta_e<=0:
                m0=m1
                Em0=Em1
                # plt.plot(m0[0],m0[1],Em0,'b*')
            else:
                r=np.random.rand(1)
                if (P > r):
                    m0=m1
                    Em0=Em1
                    # plt.plot(m0[0],m0[1],Em0,'b*') 
            muestreo_plot.append([m0[0],m0[1],Em0])
        if(T >= 2):
            T=T-1
        else:
            T=T-0.2
    S_0.append(m0[0])
    T_0.append(m0[1])

print ("S_0=",np.mean(S_0))
print ("T_0=",np.mean(T_0))
print ("Varianza S", np.var(S_0))
print ("varianza T", np.var(T_0))
print (len(muestreo_plot))
muestreo_plot_arr = np.array(muestreo_plot)
fig,ax = plt.subplots(subplot_kw = {"projection":"3d"}, figsize = (8,8))
ax.plot_trisurf (muestreo_plot_arr[:,0], muestreo_plot_arr[:,1], muestreo_plot_arr[:,2],alpha =.5)
#surf = ax.plot_surface(muestreo_plot_arr[:,0], muestreo_plot_arr[:,1],Z,cmap=plt.matplotlib.cm.get_cmap('hot'),alpha=0.5,linewidth=5)
#plt.plot(m0[0],m0[1],Em0,'mo') 
ax.plot(m0[0],m0[1],Em0,'mo',label='Punto final')
ax.legend(loc='upper left')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_zlim(0,1)
plt.show()        
