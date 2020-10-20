import numpy as np
from scipy.linalg import solve
from scipy.linalg import eigh
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
T = 300.
ni = 2.86e25*np.exp(-const.e*0.56/(const.k*T))
M1 = int(N*8./66.)
M2 = int(N*58./66.)

mzz = 0.19*const.m_e
Lx = 100e-9
Ly = 100e-9

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

def make_b(Vg,phi,n=np.zeros(N),feedback_on=False):
    b = np.zeros(N)

    b[M1+1:M2] = const.e*Nacc*h
    b[M1] = const.e*Nacc*h/2.
    b[M2] = const.e*Nacc*h/2.

    if feedback_on :
        b += const.e*n*h
    else :
        b[M1+1:M2] += 2.*const.e*ni*np.sinh(const.e*phi[M1+1:M2]/(const.k*T))*h
        b[M1] += const.e*ni*np.sinh(const.e*phi[M1]/(const.k*T))*h
        b[M2] += const.e*ni*np.sinh(const.e*phi[M2]/(const.k*T))*h

    b[0] = phi0+Vg; b[N-1] = phi0+Vg

    return b

def make_db(phi):
    db = np.zeros((N,N))

    np.fill_diagonal(db[M1+1:M2,M1+1:M2],2.*const.e*ni*np.cosh(const.e*phi[M1+1:M2]/(const.k*T))*const.e/(const.k*T)*h)
    db[M1,M1] = const.e*ni*np.cosh(const.e*phi[M1]/(const.k*T))*const.e/(const.k*T)*h
    db[M2,M2] = const.e*ni*np.cosh(const.e*phi[M2]/(const.k*T))*const.e/(const.k*T)*h

    return db

def Newton(dx_list,phi_,Vg):
    H = make_mat()

    dx = 1.; phi = phi_
    while dx > 1e-7:
        r = np.matmul(H,phi)-make_b(Vg,phi)
        J = H-make_db(phi)
        dphi = solve(J,-r)
        phi += dphi
        dx = np.linalg.norm(dphi)
        dx_list.append(dx)

    return phi

def make_mat_Ham():
    H = np.zeros((N-2,N-2))

    np.fill_diagonal(H,-2.)
    np.fill_diagonal(H[1:],1)
    np.fill_diagonal(H[:,1:],1)

    return H/h

def FD_dist(E):
    return 1./(1.+np.exp(E/(const.k*T)))

def ele_density(phi):
    gap = np.full(N-2,4.5)
    gap[M1-1:M2] = 0.56

    mzz_list = np.full(N-2,0.58*const.m_e)
    mzz_list[M1-1:M2] = mzz

    Ham = make_mat_Ham()-2.*np.diag(mzz_list)*h/const.hbar**2*const.e*np.diag(gap-phi[1:N-1])
    E,V = eigh(-Ham*const.hbar**2/(2.*mzz*h))

    V = V*np.sqrt(1./h)
    v0 = np.zeros((1,N-2))
    V = np.concatenate((v0,V),axis=0)
    V = np.concatenate((V,v0),axis=0)

    #z = np.linspace(0,L,N)
    #print(simps(V[:,0]**2,z,h))

    n_list = np.sum(V**2*FD_dist(E),axis=1)

    return 2.*n_list/(Lx*Ly), V

def self_consistency(Vg,phi_,n_qn):
    dx = 1.
    while dx>1e-7:
        phi = solve(make_mat(),make_b(Vg,phi_,n_qn,True))
        n_qn,psi = ele_density(phi)
        dphi = phi-phi_; phi_ = phi
        dx = np.linalg.norm(dphi); print(dx)
    return n_qn

phi_ = np.full(N,phi0)
z = np.linspace(0,L,N)

plt.figure(figsize=(10,10))

Vg = 0.0
dx_list=[]
phi = Newton(dx_list,phi_,Vg)
n_qn,psi = ele_density(phi)

plt.xlabel('z',fontsize=18)
plt.ylabel('$|\psi_m(z)|^2$',fontsize=18)
for i in range(4):
    plt.plot(z,psi[:,i]**2,c=plt.cm.coolwarm(i/4.),label='$m=$'+str(i+1))
plt.legend(loc='best',fontsize=18)

plt.show()

plt.figure(figsize=(10,10))

plt.xlabel('z',fontsize=18)
plt.ylabel('n(z)',fontsize=18)
plt.plot(z,n_qn,'r-',label='$V_G=$'+str(Vg))
plt.legend(loc='best',fontsize=18)

plt.show()

fig,axs = plt.subplots(1,2,figsize=(16,8))

for Vg in np.linspace(0,1,11):
    dx_list=[]
    phi = Newton(dx_list,phi_,Vg)
    n_qn,psi = ele_density(phi)
    n_cl = np.zeros(N); n_cl[M1:M2+1] =  ni*np.exp(const.e*phi[M1:M2+1]/(const.k*T))
    
    n_qn = self_consistency(Vg,phi,n_qn)

    axs[0].plot(z,n_qn,c=plt.cm.jet(Vg))
    axs[1].plot(z,n_cl,c=plt.cm.jet(Vg),label='$V_G=$'+str(format(Vg,'1.1f')))

axs[0].set_xlabel('z',fontsize=18)
axs[0].set_ylabel('n(z)',fontsize=18)
legend = axs[0].legend(loc='best',fontsize=18,frameon=False)
legend.set_title('Quantum',prop={'size':18})

axs[1].set_xlabel('z',fontsize=18)
axs[1].set_ylabel('n(z)',fontsize=18)
legend = axs[1].legend(loc='best',fontsize=18,frameon=False)
legend.set_title('Semi-classical',prop={'size':18})

plt.show()

