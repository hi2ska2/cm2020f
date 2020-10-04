import numpy as np
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 500
L = 6.6e-9
e1 = 3.9*const.epsilon_0
e2 = 11.7*const.epsilon_0
Nacc = 1e24
phi0 = 0.33374
ni = 2.86e25*np.exp(-const.e*0.56/(const.k*300))
M1 = int(N*8./66.)
M2 = int(N*58./66.)

def make_mat():
    H = np.zeros((N,N))

    np.fill_diagonal(H[:M1,:],-2.*e1)
    np.fill_diagonal(H[1:M1,:],e1)
    np.fill_diagonal(H[:M1,1:],e1)

    H[M1,M1] = -e1-e2; H[M1,M1+1] = e2; H[M1,M1-1] = e1

    np.fill_diagonal(H[M1+1:,M1+1:],-2.*e2)
    np.fill_diagonal(H[M1+1:,M1+2:],e2)
    np.fill_diagonal(H[M1+1:,M1:],e2)

    H[M2,M2] = -e2-e1; H[M2,M2+1] = e1; H[M2,M2-1] = e2

    np.fill_diagonal(H[M2+1:,M2+1:],-2.*e1)
    np.fill_diagonal(H[M2+1:,M2+2:],e1)
    np.fill_diagonal(H[M2+1:,M2:],e1)

    H[0,0] = 1.; H[N-1,N-1] = 1.
    H[0,1] = 0.; H[N-1,N-2] = 0.
    
    return H

def make_b(Vg,phi=np.zeros(N)):
    b = np.zeros(N)
    h = L/(N-1)

    b[M1+1:M2] = const.e*Nacc*h**2
    b[M1] = const.e*Nacc*h**2/2.
    b[M2] = const.e*Nacc*h**2/2.

    b += 2.*const.e*ni*np.sinh(const.e*phi/(const.k*300.))*h**2
    b[0] = phi0+Vg; b[N-1] = phi0+Vg

    return b

n = np.zeros(N)
p = np.zeros(N)
H = make_mat()
Vg = 0.1

b = make_b(Vg)
phi = solve(H,b)
n[M1+1:M2] = ni*np.exp(const.e*phi[M1+1:M2]/(const.k*300))

fig,axs = plt.subplots(1,2,figsize=(16,8))
x = np.linspace(0,L,N)
axs[0].plot(x,phi,'y-',lw=2,ms=7,mfc='None',label='initial $\phi(x)$')
axs[1].plot(x,n,'y-',lw=2,ms=7,mfc='None',label='initial $n(x)$')

b = make_b(Vg,phi)
phi = solve(H,b)
n[M1+1:M2] = ni*np.exp(const.e*phi[M1+1:M2]/(const.k*300))

axs[0].plot(x,phi,'b-.',lw=2,ms=7,mfc='None',label='updated $\phi(x)$')
axs[1].plot(x,n,'b-.',lw=2,ms=7,mfc='None',label='updated $n(x)$')

dx = 1.
while np.abs(dx) > 1e-6:
    dx = phi[N/2]
    b = make_b(Vg,phi)
    phi = solve(H,b)
    dx -= phi[N/2]

n[M1+1:M2] = ni*np.exp(const.e*phi[M1+1:M2]/(const.k*300))
axs[0].plot(x,phi,'r:',lw=2,ms=7,mfc='None',label='final $\phi(x)$')
axs[1].plot(x,n,'r:',lw=2,ms=7,mfc='None',label='final $n(x)$')

axs[0].set_xlabel('$x$',fontsize=18)
axs[0].set_ylabel(r'$\phi(x)$',fontsize=18)
legend = axs[0].legend(loc='best',frameon = False, fontsize=18)
legend.set_title('$V_G=$'+str(Vg)+'V',prop={'size':18})

axs[1].set_xlabel('$x$',fontsize=18)
axs[1].set_ylabel(r'$n(x)$',fontsize=18)

fig.tight_layout()
plt.show()

fig,axs = plt.subplots(1,1,figsize=(10,10))

for i in range(8):
    b = make_b(0.3+i*0.1)
    phi = solve(H,b)
    b = make_b(0.3+i*0.1,phi)
    phi = solve(H,b)

    axs.plot(x,phi,'C'+str(i),lw=1,ms=7,mfc='None',label='$V_G=$'+str(0.3+i*0.1)+'V')

axs.set_xlabel('$x$',fontsize=18)
axs.set_ylabel(r'$\phi(x)$',fontsize=18)
axs.legend(loc='upper left',bbox_to_anchor=(1,1),frameon = False, fontsize=18)
axs.set_yscale('symlog')

fig.tight_layout()
plt.show()

