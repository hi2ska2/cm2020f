import numpy as np
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt

L = 6.6e-9
e1 = 3.9*const.epsilon_0
e2 = 11.7*const.epsilon_0
Nacc = 1e23

def make_prob(N):
    H = np.zeros((N,N))
    b = np.zeros(N)
    h = L/(N-1)

    M1 = int(N*8./66.)
    M2 = int(N*58./66.)

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
    
    b[M1+1:M2] = const.e*Nacc*h**2
    b[M1] = const.e*Nacc*h**2/2.
    b[M2] = const.e*Nacc*h**2/2.

    return H,b

def exact(x):
    if x < 8./66.*L:
        return -const.e*Nacc/(2.*e1)*50./66.*L*x
    elif x < 58./66.*L:
        return const.e*Nacc/(2.*e2)*(x-33./66.*L)**2-const.e*Nacc/2.*(50./66.*8./66./e1+(25./66.)**2/e2)*L**2
    else:
        return const.e*Nacc/(2.*e1)*50./66.*L*(x-L)

def get_exact(x):
    return np.vectorize(exact)(x)

fig,axs = plt.subplots(1,2,figsize=(16,8))
for n,marker,color in [(50,'-','y'),(500,'--','b'),(2000,':','r')]:
    H,b = make_prob(n)
    phi = solve(H,b)

    xn = np.linspace(0,L,n)
    axs[0].plot(xn,phi,color+marker,lw=2.5,ms=7,mfc='None',label='$N=$'+str(n))
    axs[1].plot(xn,np.abs(get_exact(xn)-phi),color+marker,lw=1,ms=7,mfc='None')

axs[0].set_xlabel('$x$',fontsize=18)
axs[0].set_ylabel(r'$\phi(x)$',fontsize=18)
legend = axs[0].legend(loc='best',frameon = False, fontsize=18)
legend.set_title('$N_{acc}=$'+str(Nacc),prop={'size':18})
axs[0].set_xlim(0,L)
axs[0].set_ylim(top=0)

axs[1].set_xlabel('$x$',fontsize=18)
axs[1].set_ylabel(r'$|\phi(x)-\phi_{ext.}(x)|$',fontsize=18)
axs[1].ticklabel_format(style='sci',scilimits=(-4,4),axis='both')
axs[1].set_xlim(0,L)
axs[1].set_ylim(bottom=0)

fig.tight_layout()
plt.show()

