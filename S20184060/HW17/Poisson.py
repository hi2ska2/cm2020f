import numpy as np
import copy
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

Nx = 61
Ny = 15
N = Nx*Ny
h = 0.5e-9
Np = -1e26
Ni = 2.86e25*np.exp(-const.e*0.56/(const.k*300))
e1 = 3.9*const.epsilon_0
e2 = 11.7*const.epsilon_0
phi0 = 0.33374
phi1 = const.k*300./const.e*np.arcsinh(np.abs(Np)/Ni/2.)
My1 = 3
My2 = 12
Mx1 = 20
Mx2 = 40

class SP_Matrix:
    def __init__(self):
        self.H_row = np.array([])
        self.H_col = np.array([])
        self.H_data = np.array([])

    def __getitem__(self,rcd):
        row,col,data = rcd
        self.H_row = np.append(self.H_row,row)
        self.H_col = np.append(self.H_col,col)
        self.H_data = np.append(self.H_data,data)

    def csr_form(self):
        return (self.H_data,(self.H_row,self.H_col))

    def to_array(self):
	csr_H = csr_matrix(self.csr_form(),shape=(N,N))
	return csr_H.toarray()

class Poisson_eq_2D:
    def __init__(self,Vg):
        self.H = SP_Matrix()
        self.b = np.zeros(N)
        
	############
	## Vortex ##
	############

        self.H[0,0,-e1]; self.H[0,1,0.5*e1]; self.H[0,Nx,0.5*e1]
        self.H[Nx-1,Nx-1,-e1]; self.H[Nx-1,Nx-2,0.5*e1]; self.H[Nx-1,2*Nx-1,0.5*e1]
        self.H[(Ny-1)*Nx,(Ny-1)*Nx,-e1]; self.H[(Ny-1)*Nx,(Ny-1)*Nx+1,0.5*e1]; self.H[(Ny-1)*Nx,(Ny-2)*Nx,0.5*e1]
        self.H[N-1,N-1,-e1]; self.H[N-1,N-2,0.5*e1]; self.H[N-1,(Ny-1)*Nx-1,0.5*e1]

	#######################
	## Edge at Interface ##
	#######################

        self.H[(My1-1)*Nx,(My1-1)*Nx,-(e1+e2)]; self.H[(My1-1)*Nx,(My1-1)*Nx+1,(e1+e2)/2.]; self.H[(My1-1)*Nx,(My1-2)*Nx,0.5*e1]; self.H[(My1-1)*Nx,My1*Nx,0.5*e2]
        self.H[My1*Nx-1,My1*Nx-1,-(e1+e2)]; self.H[My1*Nx-1,My1*Nx-2,(e1+e2)/2.]; self.H[My1*Nx-1,(My1-1)*Nx-1,0.5*e1]; self.H[My1*Nx-1,(My1+1)*Nx-1,0.5*e2]
        self.H[My2*Nx,My2*Nx,-(e1+e2)]; self.H[My2*Nx,My2*Nx+1,(e1+e2)/2.]; self.H[My2*Nx,(My2-1)*Nx,0.5*e2]; self.H[My2*Nx,(My2+1)*Nx,0.5*e1]
        self.H[(My2+1)*Nx-1,(My2+1)*Nx-1,-(e1+e2)]; self.H[(My2+1)*Nx-1,(My2+1)*Nx-2,(e1+e2)/2.]; self.H[(My2+1)*Nx-1,My2*Nx-1,0.5*e2]; self.H[(My2+1)*Nx-1,(My2+2)*Nx-1,0.5*e1]

	########################
	## Bottom Oxide Layer ##
	########################

        for i in range(0,My1*Nx):
            ix = i%Nx
            iy = i/Nx
            if ix==0 and not(iy==0 or iy==My1-1) : 
                self.H[i,i,-2.*e1]; self.H[i,i+1,e1]; self.H[i,i+Nx,0.5*e1]; self.H[i,i-Nx,0.5*e1]
            elif ix==Nx-1 and not(iy==0 or iy==My1-1) : 
                self.H[i,i,-2.*e1]; self.H[i,i-1,e1]; self.H[i,i+Nx,0.5*e1]; self.H[i,i-Nx,0.5*e1]
            elif iy==0 and not(ix==0 or ix==Nx-1) :
                if ix<Mx1 or ix>Mx2 :
                    self.H[i,i,-2.*e1]; self.H[i,i+Nx,e1]; self.H[i,i-1,0.5*e1]; self.H[i,i+1,0.5*e1]
                else :
                    self.H[i,i,1.]; self.b[i] = phi0+Vg
            elif iy==My1-1 and not(ix==0 or ix==Nx-1) :
                self.H[i,i,-2.*(e1+e2)]; self.H[i,i-1,(e1+e2)/2.]; self.H[i,i+1,(e1+e2)/2.]; self.H[i,i-Nx,e1]; self.H[i,i+Nx,e2]
            elif not(ix==0 or ix==Nx-1) and not(iy==0 or iy==My1-1) :
                self.H[i,i,-4.*e1]; self.H[i,i-1,e1]; self.H[i,i+1,e1]; self.H[i,i-Nx,e1]; self.H[i,i+Nx,e1]
        
	###################
	## Silicon Layer ##
	###################

        for i in range(My1*Nx,My2*Nx):
            ix = i%Nx
            iy = i/Nx
            if ix==0 : 
                self.H[i,i,1.]; self.b[i] = phi1
            elif ix==Nx-1 : 
                self.H[i,i,1.]; self.b[i] = phi1
            else :
                self.H[i,i,-4.*e2]; self.H[i,i-1,e2]; self.H[i,i+1,e2]; self.H[i,i-Nx,e2]; self.H[i,i+Nx,e2]
                if ix<Mx1 or ix>Mx2 :
                    self.b[i] = h*h*const.e*Np
                
	#####################
	## Top Oxide Layer ##
	#####################

        for i in range(My2*Nx,N):
            ix = i%Nx
            iy = i/Nx
            if ix==0 and not(iy==My2 or iy==Ny-1) : 
                self.H[i,i,-2.*e1]; self.H[i,i+1,e1]; self.H[i,i+Nx,0.5*e1]; self.H[i,i-Nx,0.5*e1]
            elif ix==Nx-1 and not(iy==My2 or iy==Ny-1) : 
                self.H[i,i,-2.*e1]; self.H[i,i-1,e1]; self.H[i,i+Nx,0.5*e1]; self.H[i,i-Nx,0.5*e1]
            elif iy==My2 and not(ix==0 or ix==Nx-1) :
                self.H[i,i,-2.*(e1+e2)]; self.H[i,i-1,(e1+e2)/2.]; self.H[i,i+1,(e1+e2)/2.]; self.H[i,i-Nx,e2]; self.H[i,i+Nx,e1]
            elif iy==Ny-1 and not(ix==0 or ix==Nx-1) :
                if ix<Mx1 or ix>Mx2 :
                    self.H[i,i,-2.*e1]; self.H[i,i-Nx,e1]; self.H[i,i-1,0.5*e1]; self.H[i,i+1,0.5*e1]
                else :
                    self.H[i,i,1.]; self.b[i] = phi0+Vg
            elif not(ix==0 or ix==Nx-1) and not(iy==My2 or iy==Ny-1) :
                self.H[i,i,-4.*e1]; self.H[i,i-1,e1]; self.H[i,i+1,e1]; self.H[i,i-Nx,e1]; self.H[i,i+Nx,e1]

    def make_r(self,phi_):
        for i in range(My1*Nx,My2*Nx):
            ix = i%Nx
            iy = i/Nx
            if not(ix==0 or ix==Nx-1) :
                if ix<Mx1 or ix>Mx2 :
                    self.b[i] = h*h*const.e*(Np+2.*Ni*np.sinh(const.e*phi_[i]/(300.*const.k)))
                else : 
                    self.b[i] = 2.*h*h*const.e*Ni*np.sinh(const.e*phi_[i]/(300.*const.k))

        csr_H = csr_matrix(self.H.csr_form(),shape=(N,N))
        Hphi = csr_H.dot(phi_)

        return Hphi-self.b

    def make_J(self,phi_):
        J = copy.deepcopy(self.H)
        for i in range(My1*Nx,My2*Nx):
            ix = i%Nx
            iy = i/Nx
            if not(ix==0 or ix==Nx-1) :
                J[i,i,-2.*h*h*const.e*Ni*np.cosh(const.e*phi_[i]/(300.*const.k))*const.e/(300.*const.k)] 

        return J

def Newton(Poisson,phi_,dx_list):
    dx = 1.; phi = phi_
    while dx > 1e-7:
        r = Poisson.make_r(phi)
        J = Poisson.make_J(phi)
        dphi = spsolve(csr_matrix(J.csr_form(),shape=(N,N)),-r)
        #dphi = solve(J.to_array(),-r)
        phi += dphi
        dx = np.linalg.norm(dphi); print(dx)
        dx_list.append(dx)

    return phi

def Elec_density(phi):
    n = np.zeros(N)
    n[My1*Nx:My2*Nx] = Ni*np.exp(const.e*phi[My1*Nx:My2*Nx]/(300.*const.k))
    return n

