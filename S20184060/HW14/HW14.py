import numpy as np
import copy
from scipy.linalg import solve
from scipy.integrate import simps
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 800
L = 600e-9
h = L/(N-1)
e1 = 11.7*const.epsilon_0
Np = 1e23
ni = 1.075e16
Vt = const.k*300./const.e
phi0 = Vt*np.log(Np/ni)
M1 = int(N*1./2.)

def make_mat():
    H = np.zeros((N,N))

    np.fill_diagonal(H[:,:],-2.*e1/h)
    np.fill_diagonal(H[1:,:],e1/h)
    np.fill_diagonal(H[:,1:],e1/h)

    H[0,0] = 1.; H[N-1,N-1] = 1.
    H[0,1] = 0.; H[N-1,N-2] = 0.
    
    return H

def make_b(Vg,phi,n=[],p=[],feedback_on=False):
    b = np.zeros(N)
    b[:M1] = const.e*Np*h
    b[M1+1:] = -const.e*Np*h

    if feedback_on :
        b += const.e*(n-p)*h
    else :
        b += 2.*const.e*ni*np.sinh(phi/Vt)*h

    b[0] = -phi0; b[N-1] = phi0+Vg

    return b

def make_db(phi,n=[],p=[],feedback_on=False):
    db = np.zeros((N,N))

    if feedback_on :
        np.fill_diagonal(db,2.*const.e*(n+p)*np.cosh(phi/Vt)/Vt*h)
    else :
        np.fill_diagonal(db,2.*const.e*ni*np.cosh(phi/Vt)/Vt*h)

    return db

def make_Je(phi,n):
    Je = np.zeros((N,N))
    r = np.zeros(N)

    dphi = np.zeros(N); dphi[:N-1] = (phi[1:]-phi[:N-1])/h
    dn = np.zeros(N); dn[:N-1] = (n[1:]-n[:N-1])/h
    nav = np.zeros(N); nav[:N-1] = (n[1:]+n[:N-1])/2.

    Jn = nav*dphi-Vt*dn
    dJn = 0.5*dphi+Vt/h
    dJn_ = 0.5*dphi-Vt/h

    r[1:] = Jn[1:]-Jn[:N-1]
    np.fill_diagonal(Je[1:,1:],dJn[1:]-dJn_[:N-1])
    np.fill_diagonal(Je[:,1:],dJn_)
    np.fill_diagonal(Je[1:,:],-dJn)

    r[0] = n[0]-ni**2/Np; r[N-1] = n[N-1]-Np
    Je[0,1] = 0.; Je[0,0] = 1.
    Je[N-1,N-2] = 0.; Je[N-1,N-1] = 1.

    return r,Je

def make_Jh(phi,p):
    Jh = np.zeros((N,N))
    r = np.zeros(N)

    dphi = np.zeros(N); dphi[:N-1] = (phi[1:]-phi[:N-1])/h
    dp = np.zeros(N); dp[:N-1] = (p[1:]-p[:N-1])/h
    pav = np.zeros(N); pav[:N-1] = (p[1:]+p[:N-1])/2.

    Jp = -pav*dphi-Vt*dp
    dJp = -0.5*dphi+Vt/h
    dJp_ = -0.5*dphi-Vt/h

    r[1:] = Jp[1:]-Jp[:N-1]
    np.fill_diagonal(Jh[1:,1:],dJp[1:]-dJp_[:N-1])
    np.fill_diagonal(Jh[:,1:],dJp_)
    np.fill_diagonal(Jh[1:,:],-dJp)

    r[0] = p[0]-Np; r[N-1] = p[N-1]-ni**2/Np
    Jh[0,1] = 0.; Jh[0,0] = 1.
    Jh[N-1,N-2] = 0.; Jh[N-1,N-1] = 1.

    return r,Jh

def Newton(phi_,Vg):
    H = make_mat()

    dx = 1.; phi = phi_
    while dx > 1e-10:
        r = np.matmul(H,phi)-make_b(Vg,phi)
        J = H-make_db(phi)
        dphi = solve(J,-r)
        phi += dphi
        dx = np.linalg.norm(dphi)

    return phi

def Newton_feedback(phi_,Vg,n,p):
    H = make_mat()

    dx_list = []
    dx = 1.; phi = phi_
    while dx > 1e-7:
        r = np.matmul(H,phi)-make_b(Vg,phi,n,p,True)
        J = H-make_db(phi,n,p,True)
        dphi = solve(J,-r)
        phi += dphi
        dx = np.linalg.norm(dphi); #print('Nt fb : '+str(dx))
	dx_list.append(dx)

    #plt.plot(dx_list); plt.show()

    return phi

def Newton_DD(dx_list,phi_,Vg):
    H = make_mat()

    dx = 1.; phi = phi_
    n = ni*np.exp(phi/Vt)
    p = ni*np.exp(-phi/Vt)

    while dx > 1e-7:
        r,Je = make_Je(phi,n); #plt.plot(r); plt.show()
        dn = solve(Je,-r)
        n += dn

        r,Jh = make_Jh(phi,p); #plt.plot(r); plt.show()
        dp = solve(Jh,-r)
        p += dp

        phi = Newton_feedback(phi,Vg,n,p)

        dx = np.max(np.abs(dn/n)); print('Nt DD : '+str(dx)+'\n')
        dx_list.append(dx)

    return n,p

phi_ = np.full(N,phi0)
x = np.linspace(0,L,N)

phi = Newton(phi_,0.0)

plt.plot(x,phi,'b-')

plt.xlabel('x',fontsize=18)
plt.ylabel('$\phi(x)$',fontsize=18)

plt.show()

dx_list = []
n,p = Newton_DD(dx_list,phi,0.0)

plt.plot(x,p,'r-',mfc='None',label='w/ Drift-Diffusion')
plt.plot(x[::10],ni*np.exp(-phi[::10]/Vt),'bo',mfc='None',label='w/o Drift-Diffusion')

plt.xlabel('x',fontsize=18)
plt.ylabel('p(x)',fontsize=18)
plt.yscale('log')
plt.legend(loc='upper right',frameon=False,fontsize=15)

plt.show()

