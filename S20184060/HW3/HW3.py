import numpy as np
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt

L = 1.

def make_prob(N):
    h = L/(N-1)
    H = np.zeros((N,N))
    b1 = np.zeros(N)
    b2 = np.zeros(N)

    np.fill_diagonal(H,-2.)
    np.fill_diagonal(H[1:],1.)
    np.fill_diagonal(H[:,1:],1.)
    H = H/h**2

    H[0,0] = 1.; H[N-1,N-1] = 1.
    H[0,1] = 0.; H[N-1,N-2] = 0.
    b1[0] = 1.; b1[N-1] = -1. # prob 1
    b2[(N-1)/2] = 1./h # prob 2

    return H,b1,b2

def prob1(x):
    return -2.*x+1.

def prob2(x):
    return 0.5*np.abs(x-0.5)-0.25

N = 499
H,b1,b2 = make_prob(N)
phi1 = solve(H,b1)
phi2 = solve(H,b2)

x = np.linspace(0,L,N)
fig,axs = plt.subplots(1,2,figsize=(16,8))
axs[0].plot(x,phi1,'r-',lw=0.5,ms=7,mfc='None',label='Prob.1')
axs[0].plot(x,phi2,'b-',lw=0.5,ms=7,mfc='None',label='Prob.2')

for n,marker,intv in [(5,'o',1),(49,'s',5),(499,'^',10)]:
    H,b1,b2 = make_prob(n)
    phi1 = solve(H,b1)
    phi2 = solve(H,b2)

    xn = np.linspace(0,L,n)
    axs[1].plot(xn,np.abs(prob1(xn)-phi1),'r-',lw=0.5,ms=7,mfc='None')
    axs[1].plot(xn,np.abs(prob2(xn)-phi2),'b-.',lw=0.5,ms=7,mfc='None')
    axs[1].plot(xn[::intv],np.abs(prob1(xn)-phi1)[::intv],'r'+marker,lw=0.5,ms=7,label='N='+str(n)+' prob.1')
    axs[1].plot(xn[::intv],np.abs(prob2(xn)-phi2)[::intv],'b'+marker,lw=0.5,ms=7,mfc='None',label='N='+str(n)+' prob.2')

axs[0].set_xlabel('$x$',fontsize=18)
axs[0].set_ylabel(r'$\phi(x)$',fontsize=18)
legend = axs[0].legend(loc='best',frameon = False, fontsize=18)
legend.set_title('$N=$'+str(N),prop={'size':18})
axs[0].set_ylim(-1,1)
axs[0].set_xlim(0,1)

axs[1].set_xlabel('$x$',fontsize=18)
axs[1].set_ylabel(r'$|\phi(x)-\phi_{ext.}(x)|$',fontsize=18)
legend = axs[1].legend(loc='best',frameon = False, fontsize=18)
axs[1].set_ylim(0,8e-14)
axs[1].set_xlim(0,1)

fig.tight_layout()
plt.show()
