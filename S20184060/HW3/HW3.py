import numpy as np
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 499 # must be odd
L = 1.
h = L/(N-1)

H = np.zeros((N,N))
np.fill_diagonal(H,-2.)
np.fill_diagonal(H[1:],1.)
np.fill_diagonal(H[:,1:],1.)
H = H/h**2

H[0,0] = 1.; H[N-1,N-1] = 1.
H[0,1] = 0.; H[N-1,N-2] = 0.

b1 = np.zeros(N); b1[0] = 1.; b1[N-1] = -1. # prob 1
b2 = np.zeros(N); b2[(N-1)/2] = 1./h # prob 2

phi1 = solve(H,b1)
phi2 = solve(H,b2)

x = np.linspace(0,L,N)
fig,axs = plt.subplots(1,1,figsize=(10,8))
axs.plot(x,phi1,'r-',lw=0.5,ms=7,mfc='None',label='Prob.1')
axs.plot(x,phi2,'b-',lw=0.5,ms=7,mfc='None',label='Prob.2')

axs.set_xlabel('$x$',fontsize=18)
axs.set_ylabel(r'$\phi(x)$',fontsize=18)
legend = axs.legend(loc='best',frameon = False, fontsize=18)
legend.set_title('$N=$'+str(N),prop={'size':18})
axs.set_ylim(-1,1)
axs.set_xlim(0,1)

fig.tight_layout()
plt.show()
