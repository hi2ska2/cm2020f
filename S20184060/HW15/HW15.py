import numpy as np
import copy
from scipy.linalg import solve
from scipy.integrate import simps
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 601
L = 60e-9
h = L/(N-1)
e1 = 11.7
Np1 = 5e23
Np2 = 2e21
ni = 1.075e16
Vt = const.k*300./const.e
phi0 = Vt*np.log(Np1/ni)
M1 = int(N*1./6.)
M2 = int(N*5./6.)

def make_R(phi,n):
    H = np.zeros((N,N))

    np.fill_diagonal(H,-2.*e1)
    np.fill_diagonal(H[1:],e1)
    np.fill_diagonal(H[:,1:],e1)

    Np_list = np.full(N,Np1); Np_list[M1:M2+1] = Np2
    Rphi = np.matmul(H,phi)+const.e/const.epsilon_0*(Np_list-n)*h**2
    Rphi[0] = phi[0]-phi0; Rphi[N-1] = phi[N-1]-phi0
    
    dphi = np.zeros(N); dphi[:N-1] = (phi[1:]-phi[:N-1])/h
    dn = np.zeros(N); dn[:N-1] = (n[1:]-n[:N-1])/h
    nav = np.zeros(N); nav[:N-1] = (n[1:]+n[:N-1])/2.

    Jn = nav*dphi-Vt*dn
    Rn = np.zeros(N); Rn[1:] = Jn[1:]-Jn[:N-1]
    Rn[0] = n[0]-Np1; Rn[N-1] = n[N-1]-Np1

    R = np.zeros(2*N)
    R[0::2] = Rphi
    R[1::2] = Rn

    return R

def make_J(phi,n):
    dphi_Rphi = np.zeros((N,N))
    dphi_Rn = np.zeros((N,N))
    dn_Rphi = np.zeros((N,N))
    dn_Rn = np.zeros((N,N))

    np.fill_diagonal(dphi_Rphi,-2.*e1)
    np.fill_diagonal(dphi_Rphi[1:],e1)
    np.fill_diagonal(dphi_Rphi[:,1:],e1)
    dphi_Rphi[0,1] = 0.; dphi_Rphi[0,0] = 1.
    dphi_Rphi[N-1,N-2] = 0.; dphi_Rphi[N-1,N-1] = 1.

    np.fill_diagonal(dn_Rphi,-const.e/const.epsilon_0*h**2)
    dn_Rphi[0,1] = 0.; dn_Rphi[0,0] = 0.
    dn_Rphi[N-1,N-2] = 0.; dn_Rphi[N-1,N-1] = 0.

    dphi = np.zeros(N); dphi[:N-1] = (phi[1:]-phi[:N-1])/h
    nav = np.zeros(N); nav[:N-1] = (n[1:]+n[:N-1])/2.

    np.fill_diagonal(dphi_Rn[1:,1:],-(nav[1:]+nav[:N-1])/h)
    np.fill_diagonal(dphi_Rn[1:],nav/h)
    np.fill_diagonal(dphi_Rn[:,1:],nav/h)
    dphi_Rn[0,1] = 0.; dphi_Rn[0,0] = 0.
    dphi_Rn[N-1,N-2] = 0.; dphi_Rn[N-1,N-1] = 0.

    dJn = 0.5*dphi+Vt/h
    dJn_ = 0.5*dphi-Vt/h

    np.fill_diagonal(dn_Rn[1:,1:],dJn[1:]-dJn_[:N-1])
    np.fill_diagonal(dn_Rn[:,1:],dJn_)
    np.fill_diagonal(dn_Rn[1:,:],-dJn)
    dn_Rn[0,1] = 0.; dn_Rn[0,0] = 1.
    dn_Rn[N-1,N-2] = 0.; dn_Rn[N-1,N-1] = 1.

    J = np.zeros((2*N,2*N))
    J[0::2,0::2] = dphi_Rphi
    J[0::2,1::2] = dn_Rphi
    J[1::2,0::2] = dphi_Rn
    J[1::2,1::2] = dn_Rn

    return J

def make_mat():
    H = np.zeros((N,N))

    np.fill_diagonal(H[:,:],-2.*e1/h)
    np.fill_diagonal(H[1:,:],e1/h)
    np.fill_diagonal(H[:,1:],e1/h)

    H[0,0] = 1.; H[N-1,N-1] = 1.
    H[0,1] = 0.; H[N-1,N-2] = 0.
    
    return H

def make_b(phi):
    b = np.full(N,-Np1*const.e/const.epsilon_0*h)
    b[M1:M2+1] = -Np2*const.e/const.epsilon_0*h

    b += 2.*const.e/const.epsilon_0*ni*np.sinh(phi/Vt)*h

    b[0] = phi0; b[N-1] = phi0

    return b

def make_db(phi):
    db = np.zeros((N,N))

    np.fill_diagonal(db,2.*const.e/const.epsilon_0*ni*np.cosh(phi/Vt)/Vt*h)

    db[0] = 0.; db[N-1] = 0.

    return db

def Newton(phin):
    dx=1.
    while dx > 1e-10:
        R = make_R(phin[0::2],phin[1::2]); #plt.plot(R[1::2]); plt.show()
        J = make_J(phin[0::2],phin[1::2]); #plt.imshow(J); plt.show()
        dphin = solve(J,-R); #plt.plot(dphin[1::2]); plt.show()
        phin += dphin
        dx = np.max(np.abs(dphin/phin)); print(dx)

    return phin

def Newton_NonLin(phi):
    H = make_mat()

    dx = 1.
    while dx > 1e-10:
        r = np.matmul(H,phi)-make_b(phi)
        J = H-make_db(phi)
        dphi = solve(J,-r)
        phi += dphi
        dx = np.linalg.norm(dphi); print(dx)

    return phi

phin_ = np.full(2*N,phi0); phin_[1::2] = ni*np.exp(phi0/Vt)
phi_ = np.full(N,phi0)
x = np.linspace(0,L,N)

phin = Newton(phin_)
phi = Newton_NonLin(phi_)

fig,axs = plt.subplots(1,2,figsize=(15,8))

axs[0].plot(x,phin[0::2],'r-',label='Poisson+Continuity')
axs[0].plot(x[::1],phi[::1],'bo',mfc='None',label='Non-linear Poisson')

axs[0].set_xlabel('x',fontsize=18)
axs[0].set_ylabel('$\phi(x)$',fontsize=18)
axs[0].legend(loc='best',fontsize=15,frameon=False)

axs[1].plot(x,phin[1::2],'r-')
axs[1].plot(x[::1],ni*np.exp(phi[::1]/Vt),'bo',mfc='None')
axs[1].set_yscale('log')

axs[1].set_xlabel('x',fontsize=18)
axs[1].set_ylabel('n(x)',fontsize=18)

plt.show()

