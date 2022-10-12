from math import pi
from math import cos
from math import sqrt
from math import exp
import numpy as np
import warnings
warnings.filterwarnings('ignore')
import scipy.io as spy

#import matplotlib.pyplot as plt

def matyas(xx):
    x1 = xx[0]
    x2 = xx[1]
    term1 = 0.26 * (x1**2 + x2**2);
    term2 = -0.48*x1*x2;
    y = term1 + term2;
    return y

def ackley(xx, a=20, b=0.2, c=2*pi):
    d=len(xx)
    sum1 = 0.0;
    sum2 = 0.0;
    for ii in range (0,d):
        xi = xx[ii]
        sum1 = sum1 + xi**2
        sum2 = sum2 + cos(c*xi)

    term1 = -a * exp(-b*sqrt(sum1/d));
    term2 = -exp(sum2/d);
    y = term1 + term2 + a + exp(1);
    return y


def boha1(xx):
    x1 = xx[0]
    x2 = xx[1]
    term1 = x1**2
    term2 = 2*x2**2
    term3 = -0.3 * cos(3*pi*x1)
    term4 = -0.4 * cos(4*pi*x2)
    y = term1 + term2 + term3 + term4 + 0.7
    return y

def drop(xx):
    x1 = xx[0]
    x2 = xx[1]
    frac1 = 1 + cos(12*sqrt(x1**2+x2**2));
    frac2 = 0.5*(x1**2+x2**2) + 2;
    y = -frac1/frac2;
    return y

def punto1(xx):
    S = xx[0]
    T = xx[1]
    #Valores medidos
    hh=np.array([0.72,0.49,0.30,0.20,0.16,0.12])
    t=np.array([5.0,10,20,30,40,50])
    #sigma=0.01
    #Parametros medidos
    Q=50
    d=60
    #Ecuacion nolineal
    n=len(hh)
    dif=0
    h=np.zeros([n])
    for i in range (n):
        h[i]=(Q/(4*pi*T*t[i]))*np.exp((-(d**2)*S)/(4*T*t[i]))
        dif=dif+((h[i]-hh[i])**2)
    dif=np.sqrt(dif)
    return dif

def punto2(xx):
    alpha=xx[0]
    n=xx[1]
    h=np.array([-2,-6,-8,-12,-14,-17,-26,-35,-46,-59,-71,-85,-97,-107,-118,-125,-144,-167,-182,-209,-230,-266,-321,-388,-457,-514,-599,-647])
    theta_hobs=np.array([0.440000000000000,0.440000000000000,0.440000000000000,0.440000000000000,0.440000000000000,0.440000000000000,0.420000000000000,0.420000000000000,0.410000000000000,0.370000000000000,0.350000000000000,0.330000000000000,0.320000000000000,0.300000000000000,0.290000000000000,0.290000000000000,0.280000000000000,0.260000000000000,0.250000000000000,0.240000000000000,0.220000000000000,0.210000000000000,0.200000000000000,0.180000000000000,0.170000000000000,0.160000000000000,0.150000000000000,0.140000000000000])
    theta_r= 0.09
    theta_s =0.44
    #sigma=0.02
    #Ecuacion no lineal
    m=len(h)
    dif=0
    theta_hmod=np.zeros([m])
    for i in range (m):
        theta_hmod[i]= theta_r+(theta_s - theta_r)/(1+(-alpha*h[i])**n)**(1-1/n)
        dif+=np.sqrt((theta_hmod[i]-theta_hobs[i])**2)
        #dif+=(theta_hmod[i]-theta_hobs[i])
    return dif


# import scipy.io as spy
# h= spy.loadmat('jitorres_vgdata.mat')['h']
# theta = spy.loadmat('jitorres_vgdata.mat')['theta']
# def punto2(xx):
#     alpha=xx[0]
#     n=xx[1]
#     theta_r= 0.09
#     theta_s =0.44
#     sigma=0.02
#     #Ecuacion no lineal
#     m=len(h)
#     residuales = []
#     for i in range (m):
#         h1 = h[i]
#         theta_hmod= theta_r+(theta_s - theta_r)/(1+(-alpha*h[i])**n)**(1-1/n)
#         r = theta_hmod - theta[i]
#         residuales.append(r)
#     E = np.linalg.norm(residuales,2)
#     return E