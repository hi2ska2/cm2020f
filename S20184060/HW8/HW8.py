import numpy as np
from scipy.linalg import solve
from scipy.integrate import simps
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 800
L = 6.6e-9
h = L/(N-1)
e1 = 3.9*const.epsilon_0
e2 = 11.7*const.epsilon_0
Nacc = 1e24
phi0 = 0.33374
ni = 2.86e25*np.exp(-const.e*0.56/(const.k*300))
M1 = int(N*8./66.)
M2 = int(N*58./66.)

def make_mat():
    H = np.zeros((N,N))

    np.fill_diagonal(H[:M1,:],-2.*e1/h)
    np.fill_diagonal(H[1:M1,:],e1/h)
    np.fill_diagonal(H[:M1,1:],e1/h)

    H[M1,M1] = (-e1-e2)/h; H[M1,M1+1] = e2/h; H[M1,M1-1] = e1/h

    np.fill_diagonal(H[M1+1:,M1+1:],-2.*e2/h)
    np.fill_diagonal(H[M1+1:,M1+2:],e2/h)
    np.fill_diagonal(H[M1+1:,M1:],e2/h)

    H[M2,M2] = (-e2-e1)/h; H[M2,M2+1] = e1/h; H[M2,M2-1] = e2/h

    np.fill_diagonal(H[M2+1:,M2+1:],-2.*e1/h)
    np.fill_diagonal(H[M2+1:,M2+2:],e1/h)
    np.fill_diagonal(H[M2+1:,M2:],e1/h)

    H[0,0] = 1.; H[N-1,N-1] = 1.
    H[0,1] = 0.; H[N-1,N-2] = 0.
    
    return H

def make_b(Vg,phi):
    b = np.zeros(N)

    b[M1+1:M2] = const.e*Nacc*h
    b[M1] = const.e*Nacc*h/2.
    b[M2] = const.e*Nacc*h/2.

    b[M1+1:M2] += 2.*const.e*ni*np.sinh(const.e*phi[M1+1:M2]/(const.k*300.))*h
    b[M1] += const.e*ni*np.sinh(const.e*phi[M1]/(const.k*300.))*h
    b[M2] += const.e*ni*np.sinh(const.e*phi[M2]/(const.k*300.))*h

    b[0] = phi0+Vg; b[N-1] = phi0+Vg

    return b

def make_db(phi):
    db = np.zeros((N,N))

    np.fill_diagonal(db[M1+1:M2,M1+1:M2],2.*const.e*ni*np.cosh(const.e*phi[M1+1:M2]/(const.k*300.))*const.e/(const.k*300.)*h)
    db[M1,M1] = const.e*ni*np.cosh(const.e*phi[M1]/(const.k*300.))*const.e/(const.k*300.)*h
    db[M2,M2] = const.e*ni*np.cosh(const.e*phi[M2]/(const.k*300.))*const.e/(const.k*300.)*h

    return db

def Newton(phi_list,phi_,Vg):
    H = make_mat()

    dx = 1.; phi = phi_
    while dx > 1e-7:
        r = np.matmul(H,phi)-make_b(Vg,phi)
        J = H-make_db(phi)
        dphi = solve(J,-r)
        phi += dphi
        phi_list.append(phi)
        dx = np.linalg.norm(dphi)

    return phi

n = np.zeros(N)
phi_ = np.full(N,phi0)
x = np.linspace(0,L,N)

fig,axs = plt.subplots(1,2,figsize=(16,8))

for i in range(6):
    phi_list = []
    phi = Newton(phi_list,phi_,0.2*i)
    n[M1+1:M2] = ni*np.exp(const.e*phi[M1+1:M2]/(const.k*300))

    axs[0].plot(x,phi,'C'+str(i),lw=1,mfc='None',label='$V_G=$'+str(i*0.2)+'V')
    axs[1].plot(x,n,'C'+str(i),lw=1,mfc='None')

axs[0].set_xlabel('$x$',fontsize=18)
axs[0].set_ylabel(r'$\phi(x)$',fontsize=18)
axs[0].legend(loc='center',frameon = False, fontsize=18)

axs[1].set_xlabel('$x$',fontsize=18)
axs[1].set_ylabel(r'$n(x)$',fontsize=18)
axs[1].set_ylim(bottom=0)
axs[1].set_yscale('symlog')

fig.tight_layout()
plt.show()

Integrated_n = []
Vg_list = np.linspace(0,1,100)
for v in Vg_list:
    phi_list = []
    phi = Newton(phi_list,phi_,v)
    n[M1+1:M2] = ni*np.exp(const.e*phi[M1+1:M2]/(const.k*300))
    Integrated_n.append(simps(n,x,h))

fig,axs = plt.subplots(1,1,figsize=(8,8))
axs.plot(Vg_list,Integrated_n,'r-',mfc='None',ms=7)

axs.set_xlabel('$V_G$',fontsize=18)
axs.set_ylabel('$\int n(x)\;dx$',fontsize=18)
axs.set_xlim(0,1)
#axs.set_ylim(bottom=0)
axs.set_yscale('log')

fig.tight_layout()
plt.show()

