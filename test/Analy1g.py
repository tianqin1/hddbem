import numpy as np
from Material import Material
from sympy.solvers import solve
from sympy import *
import math
import matplotlib.pyplot as plt
X = 30

fuel = Material()
fuel.set([[1.36, 0.0181, 0.0279, 1.0]])
refl = Material()
#refl.set([[0.55, 0.0181, 0.0, 1.0]])
refl.set([[0.55, 0.0127, 0.0, 1.0]])
#  B2 =-(fuel.siga(0)-fuel.nusigf(0)/keff)/fuel.dif(0)

def func_1(keff):
    B2 =-(fuel.siga(0)-fuel.nusigf(0)/keff)/fuel.dif(0)
    b=np.sqrt(abs(B2))
    return np.tan(b*X)

def func_2(keff):
    B2 =-(fuel.siga(0)-fuel.nusigf(0)/keff)/fuel.dif(0)
    b=np.sqrt(abs(B2))
    K2 = (refl.siga(0))/refl.dif(0)
    k=np.sqrt(abs(K2))
    return k/b

# x=np.arange(-1,1,0.01)
# plt.plot(x,func_1(x))
# plt.plot(x,func_2(x))
# plt.show()

# ## validation for numerical result
#B2 = 0.00186101
#keff = 1.3523357962138103
def root_finding(keff):
    B2 =-(fuel.siga(0)-fuel.nusigf(0)/keff)/fuel.dif(0)
    B=np.sqrt(abs(B2))
    K2 = (refl.siga(0))/refl.dif(0)
    K=np.sqrt(abs(K2))
    aa = np.array([[1,-np.cos(B*X),(1/B)*np.sin(B*X)],[-np.cos(B*X),1,0],[0,1,1/K]])
    return np.linalg.det(aa)
def eig():
    for i in range(100000000):
        keff=float(i/100000000)+1.35
        val =root_finding(keff)
        if(val<1.0e-8):
            print(keff," ",val)
            break
    return keff
def get_1g_flx(keff):
    bem_flx=np.loadtxt(".flux.txt")
    c1 = bem_flx[0,1]
    nnodal=1000
    width=30
    dx=width/nnodal
    x=np.zeros(nnodal)
    flux=np.zeros(nnodal)
    flux_tmp=np.zeros(nnodal)
    B2 =-(fuel.siga(0)-fuel.nusigf(0)/keff)/fuel.dif(0)
    b=np.sqrt(abs(B2))
    K2 = (refl.siga(0))/refl.dif(0)
    k=np.sqrt(abs(K2))
    print(B2)
    print(K2)
    for i in range(nnodal):
        x[i]=i/nnodal*2*30+dx
        #print(x[i])
        if(i<500):
            flux_tmp[i]=c1*np.cos(b*x[i])
        else:
            #flux[i]=0
            flux_tmp[i]=c1*np.cos(b*width)*np.exp(-k*(x[i]-width))
    c2=bem_flx[500,1]/flux_tmp[500]
    flux[0:500]=flux_tmp[0:500]
    flux[500:nnodal]=c2*flux_tmp[500:nnodal]
    plt.plot(x,flux,'-',label='Analytical')
    plt.plot(bem_flx[:,0],bem_flx[:,1],'--',label='HDD-BEM')
    plt.grid(True)
    plt.legend()
    plt.show()
root_finding(1.411019863016964)
get_1g_flx(keff=1.411019863016964)
