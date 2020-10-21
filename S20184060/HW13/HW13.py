import numpy as np
from scipy.linalg import solve
from scipy.linalg import eigh
from scipy.integrate import simps
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 300
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

mxx = 0.91*const.m_e
myy = 0.91*const.m_e
mzz = 0.19*const.m_e
Lx = 100e-9
Ly = 100e-9

#####################
#### Poisson eq. ####
#####################

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

#########################
#### Schrodinger eq. ####
#########################

def make_mat_Ham():
    H = np.zeros((N-2,N-2))

    np.fill_diagonal(H,-2.)
    np.fill_diagonal(H[1:],1)
    np.fill_diagonal(H[:,1:],1)

    return H/h

def make_db_Ham(phi):
    phi_mat = np.full((N,N),phi)
    dphi_mat = np.diag(np.full(N,0.001))
    pdp_mat = phi_mat+dphi_mat

    n_mat = np.full((N,N),ele_density(phi))
    dn_mat = np.array([ ele_density(pdp_mat[i,:]) for i in range(N) ])

    db = const.e*h*(dn_mat-n_mat).T/np.diag(dphi_mat)

    return db

def FD_dist(E):
    return Lx*Ly/(2.*np.pi*const.hbar**2)*np.sqrt(mxx*myy)*const.k*T*np.log(1.+np.exp(-E/(const.k*T)))

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

    return 2.*n_list/(Lx*Ly)

################
#### Solver ####
################

def self_consistency(dx_list,phi_,n_qn,Vg):
    dx = 1.
    while dx>1e-7:
        phi = solve(make_mat(),make_b(Vg,phi_,n_qn,True))
        n_qn = ele_density(phi)
        dphi = phi-phi_; phi_ = phi
        dx = np.linalg.norm(dphi); print(dx)
        dx_list.append(dx)

    return n_qn

def Newton_Ham(dx_list,phi_,n_qn,Vg):
    H = make_mat()
    dx=1.; phi = phi_
    while dx>1e-7:
        r = np.matmul(H,phi)-make_b(Vg,phi,n_qn,True)
        J = H-make_db_Ham(phi)
        dphi = solve(J,-r)
        phi += dphi
        n_qn = ele_density(phi)
        dx = np.linalg.norm(dphi); print(dx)
        dx_list.append(dx)

    return n_qn

phi_ = np.full(N,phi0)
z = np.linspace(0,L,N)

Vg = 0.0
dx_list=[]
phi = Newton(dx_list,phi_,Vg)
n_qn = ele_density(phi)

plt.figure(figsize=(10,10))

plt.xlabel('z',fontsize=18)
plt.ylabel('n(z)',fontsize=18)
plt.plot(z,n_qn,'r-',label='$V_G=$'+str(Vg))
plt.legend(loc='best',fontsize=18)

plt.show()

fig,axs = plt.subplots(1,2,figsize=(16,8))

for Vg in np.linspace(0,1,11):
    print("         Vg="+str(Vg)+"         ")
    
    dx_list=[]
    phi = Newton(dx_list,phi_,Vg)
    n_qn = ele_density(phi)
    #n_cl = np.zeros(N); n_cl[M1:M2+1] =  ni*np.exp(const.e*phi[M1:M2+1]/(const.k*T)) 

    dx_list1 = []
    n_qn1 = Newton_Ham(dx_list1,phi,n_qn,Vg)

    axs[0].plot(z,n_qn1,c=plt.cm.jet(Vg),label='$V_G=$'+str(format(Vg,'1.2f')))
    #axs[1].plot(z,n_cl,c=plt.cm.coolwarm(Vg),label='$V_G=$'+str(format(Vg,'1.2f')))
    axs[1].plot(dx_list1,'.-',ms=12,lw=0.5,mfc='None',c=plt.cm.jet(Vg))

axs[0].set_xlabel('z',fontsize=18)
axs[0].set_ylabel('n(z)',fontsize=18)
legend = axs[0].legend(loc='upper right',ncol=1,fontsize=14,frameon=False)
legend.set_title('Quantum',prop={'size':15})
'''
axs[1].set_xlabel('z',fontsize=18)
axs[1].set_ylabel('n(z)',fontsize=18)
legend = axs[1].legend(loc='best',fontsize=18,frameon=False)
legend.set_title('Semi-classical',prop={'size':18})
'''
axs[1].set_xlabel('iteration',fontsize=18)
axs[1].set_ylabel('$\Vert \delta\phi \Vert_\mathrm{F}$',fontsize=18)
axs[1].set_yscale('log')

fig.tight_layout()
plt.show()

fig,axs = plt.subplots(1,3,figsize=(20,6.5))

for Vg in np.linspace(0,0.24,13):
    print("         Vg="+str(Vg)+"         ")
    
    dx_list=[]
    phi = Newton(dx_list,phi_,Vg)
    n_qn = ele_density(phi)

    dx_list1 = []; dx_list2 = []
    n_qn1 = Newton_Ham(dx_list1,phi,n_qn,Vg)
    n_qn2 = self_consistency(dx_list2,phi,n_qn,Vg)

    axs[0].plot(z,n_qn1,c=plt.cm.jet(Vg/0.24),label='$V_G=$'+str(format(Vg,'1.2f')))
    axs[0].plot(z,n_qn2,'k--')
    axs[1].plot(dx_list1,'.-',ms=12,lw=0.5,mfc='None',c=plt.cm.jet(Vg/0.24))
    axs[2].plot(dx_list2,'.-',ms=12,lw=0.5,mfc='None',c=plt.cm.jet(Vg/0.24))

axs[0].set_xlabel('z',fontsize=18)
axs[0].set_ylabel('n(z)',fontsize=18)
legend = axs[0].legend(loc='upper right',ncol=2,fontsize=12,frameon=False)
legend.set_title('Quantum',prop={'size':15})

axs[1].set_xlabel('iteration',fontsize=18)
axs[1].set_ylabel('$\Vert \delta\phi \Vert_\mathrm{F}$',fontsize=18)
axs[1].set_yscale('log')

axs[2].set_xlabel('iteration',fontsize=18)
axs[2].set_ylabel('$\Vert \delta\phi \Vert_\mathrm{F}$',fontsize=18)
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].set_xlim(left=0)

fig.tight_layout()
plt.show()

