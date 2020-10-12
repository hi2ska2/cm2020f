import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt

Lz = 5e-9
mxx = 0.19*const.m_e
myy = 0.19*const.m_e
mzz = 0.91*const.m_e
T = 300.

n = np.arange(1,101)
Ez = (const.hbar*np.pi*n)**2/(2.*mzz*Lz*Lz)
Ef = np.linspace(-0.3,0.1,100).reshape((100,1))
ez = (const.e*Ef-Ez)/(const.k*T)

N_list = np.sum(1./np.pi*np.sqrt(mxx*myy)/const.hbar**2*const.k*T*np.log(1+np.exp(ez)),axis=1)

plt.figure(figsize=(10,10))
plt.plot(Ef.flatten(),N_list,'r-')

plt.xlabel('$E_F$',fontsize=18)
plt.ylabel('$n_{tot}(E_F)$',fontsize=18)
plt.yscale('log')

plt.show()
