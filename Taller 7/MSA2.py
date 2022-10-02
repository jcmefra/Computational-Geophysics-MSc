#from math import sqrt
#import heuristics
from heuristics import punto2
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

#xx=[0.0126, 1.8541]
#xx=np.array(xx)
#punto2(xx)
w = 0
for w in range(0,25):
    w+=1
    intervalo1a=0.001
    intervalo1b=0.02  # 32 ackl, 100 boha1,   5 drop, 10 matyas, 1 punto1, 10 punto2
    paso1=0.001       # 1 ackl,   1 boha1, 0.5 drop,  1 matyas, 0.1 punto1, 0.5 punto2
    intervalo2a=1
    intervalo2b=10
    paso2=1

    X,Y = np.meshgrid(np.arange(intervalo1a,intervalo1b,paso1), np.arange(intervalo2a,intervalo2b,paso2))
    k1,k2=X.shape
    Z=np.zeros((k1,k2))

    for i in range(0, k1):
        for j in range(0, k2):
            a= np.array([X[i, j], Y[i,j]])
            Z[i][j]= punto2(a)

    T=10; 
    temp1=(np.random.rand(1,1)*0.019 +0.001) #[0.01-0.02]
    temp2=(np.random.rand(1,1)*9 +1) #[1-10] 
    m0=np.array([temp1, temp2])
    m0=m0.ravel()
    Em0=punto2(m0)

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(10,10))
    #view(-17,18)#ackley
    #view(-69,67)#boha1
    #view(-161,29)#drop
    #ax.view_init(30,angle)#matya

    plt.plot(m0[0],m0[1],Em0,'ro',label='Punto inicial')

    while T > 0: # loop over the temperature T
        for i in range (100): # loop over a number of random moves/temperature
            temp1=(np.random.rand(1,1)*0.019 +0.001) 
            temp2=(np.random.rand(1,1)*9 +1)
            m1=np.array([temp1, temp2])
            m1=m1.ravel()
            Em1=punto2(m1)
            Delta_e=Em1-Em0
            P=np.exp(-Delta_e/T)
            if (Delta_e<=0):
                m0=m1
                Em0=Em1
                plt.plot(m0[0],m0[1],Em0,'b*')
            else:
                r=np.random.rand(1)
                if (P > r):
                    m0=m1
                    Em0=Em1
                    plt.plot(m0[0],m0[1],Em0,'b*') 
        if(T >= 2):
            T=T-1
        else:
            T=T-0.2
            
    print(m0)     
    break   
    # surf = ax.plot_surface(X,Y,Z,cmap=cm.coolwarm,alpha=0.3,linewidth=5)
    # ax.plot(m0[0],m0[1],Em0,'mo',label='Punto final')
    # ax.legend(loc='upper left')
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_zlabel('Z')
    # plt.show()        
