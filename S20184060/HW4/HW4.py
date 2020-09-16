import numpy as np
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt

L = 2.4e-9
e1 = 3.9
e2 = 22.

def make_prob(N):
    H = np.zeros((N,N))
    b = np.zeros(N)

    M = int(N*5./24.)

    np.fill_diagonal(H[:M,:],-2.*e1)
    np.fill_diagonal(H[1:M,:],e1)
    np.fill_diagonal(H[:M,1:],e1)

    np.fill_diagonal(H[M+1:,M+1:],-2.*e2)
    np.fill_diagonal(H[M+1:,M+2:],e2)
    np.fill_diagonal(H[M+1:,M:],e2)

    H[M,M] = -e1-e2; H[M,M+1] = e2; H[M,M-1] = e1

    H[0,0] = 1.; H[N-1,N-1] = 1.
    H[0,1] = 0.; H[N-1,N-2] = 0.
    b[N-1] = 1.

    return H,b

def exact(x):
    return np.vectorize(lambda y: y/(5./24.+19.*e1/(24.*e2))/L if y < 5./24.*L else y/(5.*e2/(24.*e1)+19./24.)/L+(1./e1-1./e2)/(1./e1+19./(5.*e2)) )(x)

print("\nExact capacitance : "+str(const.epsilon_0/(5./(24.*e1)+19./(24.*e2))/L))

fig,axs = plt.subplots(1,2,figsize=(16,8))
for n,marker,color in [(50,'-','y'),(500,'--','b'),(10000,':','r')]:
    H,b = make_prob(n)
    phi = solve(H,b)

    M = int(n*5./24.)
    print("Numerical capacitance with N="+str(n)+" : "+str(e1*const.epsilon_0*24.*phi[M]/(5.*L)))

    xn = np.linspace(0,L,n)
    axs[0].plot(xn,phi,color+marker,lw=2.5,ms=7,mfc='None',label='$N=$'+str(n))
    axs[1].plot(xn,np.abs(exact(xn)-phi),color+marker,lw=1,ms=7,mfc='None')
print(" ")

axs[0].set_xlabel('$x$',fontsize=18)
axs[0].set_ylabel(r'$\phi(x)$',fontsize=18)
axs[0].legend(loc='best',frameon = False, fontsize=18)
axs[0].set_xlim(0,L)
axs[0].set_ylim(0,1)

axs[1].set_xlabel('$x$',fontsize=18)
axs[1].set_ylabel(r'$|\phi(x)-\phi_{ext.}(x)|$',fontsize=18)
axs[1].set_xlim(0,L)
axs[1].set_ylim(0,0.006)

fig.tight_layout()
plt.show()

