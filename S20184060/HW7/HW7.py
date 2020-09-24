import numpy as np
from scipy.linalg import solve
import scipy.constants as const
import matplotlib.pyplot as plt

Np = -1e24
Ni = 1e16
phi0 = -0.5

def get_r(phi,N=Np):
    return ( N-2.*Ni*np.sinh(const.e*phi/(const.k*300.)),-2.*Ni*np.cosh(const.e*phi/(const.k*300.))*const.e/(const.k*300.) )

def get_dphi(phi,N=Np):
    return N/(2.*Ni*np.cosh(const.e*phi/(const.k*300.))*const.e/(const.k*300.))-np.tanh(const.e*phi/(const.k*300.))*const.k*300./const.e

def Newton(dphi_list,phi0,N=Np):
    dphi = 1.; phi = phi0
    while np.abs(dphi) > 1e-7:
        dphi = get_dphi(phi,N)
        phi += dphi
        dphi_list.append(dphi)

    return phi

def Exact(N=Np):
    return const.k*300./const.e*np.arcsinh(N/Ni/2.)

print('Exact phi : '+str(Exact()))

dphi_list = []
phi = Newton(dphi_list,phi0)

print('phi : '+str(phi))

plt.plot(dphi_list,'ro',ms=10,mfc='None',label='$N^+=$'+str(Np))
plt.grid()
plt.xlabel('iteration',fontsize=18)
plt.ylabel(r'$\delta\phi$',fontsize=18)
plt.legend(loc='best',frameon=False,fontsize=15)

plt.show()

phi_list = []
Np_list = [-1e24,-1e23,-1e22,-1e20,-1e18,-1e16,1e16,1e18,1e20,1e22,1e23,1e24]
for i,n in enumerate(Np_list):
    phi0 = -0.4+i*0.4/5.
    phi_list.append(Newton(dphi_list,phi0,n))

plt.plot(Np_list,phi_list,'ro',mfc='None',ms=10,label='Numerical $\phi$')
plt.plot(Np_list,Exact(np.array(Np_list)),'rx',ms=10,label='Exact $\phi$')
plt.grid()
plt.xlabel('$n_{int}$',fontsize=18)
plt.ylabel(r'$\phi$',fontsize=18)
plt.legend(loc='best',frameon=False,fontsize=15)

plt.show()

