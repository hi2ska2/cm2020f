import numpy as np
import copy
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from Poisson import *

Nx = 61
Ny = 15
N = Nx*Ny
h = 0.5e-9
T = 300.
Vt = const.k*T/const.e
Dn = 0.15*Vt
Np = -1e26
Ni = 2.86e25*np.exp(-0.56/Vt)
e1 = 3.9
e2 = 11.7
phi0 = 0.33374
phi1 = Vt*np.arcsinh(np.abs(Np)/Ni/2.)
My1 = 3
My2 = 12
Mx1 = 20
Mx2 = 40

def Bn(x):
    x_abs = np.abs(x)
    if x_abs < 2.5020000e-2 :
            sx = x*x
            return (1.-x/2.+sx/12.*(1.-sx/60.*(1.-sx/42.)))
    elif x_abs < 1.5000000e-1 :
            sx = x*x
            return (1.-x/2.+sx/12.*(1.-sx/60.*(1.-sx/42.*(1.-sx/40.*(1.-0.025252525252525252525*sx)))))
    elif x_abs > 150.01 :
            inv_exp = np.exp(-x)
            return x*inv_exp
    else :
            inv_exp = 1./(np.exp(x)-1.)
            return x/(np.exp(x)-1.)

def dBn(x):
    x_abs = np.abs(x)
    if x_abs < 2.5020000e-2 :
            sx = x*x
            return (-0.5+x/6.*(1.-sx/30.*(1.-sx/28.)))
    elif x_abs < 1.5000000e-1 :
            sx = x*x
            return (-0.5-x/6.*(1.-sx/30.*(1.-sx/28.*(1.-sx/30.*(1.-0.031565656565656565657*sx)))))
    elif x_abs > 150.01 :
            inv_exp = np.exp(-x)
            return inv_exp-x*inv_exp
    else :
            inv_exp = 1./(np.exp(x)-1.)
            return inv_exp-x*inv_exp*(inv_exp+1.)

class Drift_Diffusion_2D:
    def __init__(self,phi,n):
        self.Rp = np.zeros(N)
        self.Rn = np.zeros(N)

        self.dpRp = np.zeros((N,N))
        self.dpRn = np.zeros((N,N))
        self.dnRp = np.zeros((N,N))
        self.dnRn = np.zeros((N,N))
   
        self.R = np.zeros(2*N)
        self.J = np.zeros((2*N,2*N))

        self.phi = copy.deepcopy(phi)
        self.n = copy.deepcopy(n)

        Cvec = np.zeros(2*N)
        Cvec[0::2] = Vt; Cvec[1::2] = Np
        self.Cmat = np.diag(Cvec)

    def Jn(self,ind1,ind2):
        return self.n[ind1]*Bn((self.phi[ind1]-self.phi[ind2])/Vt)-self.n[ind2]*Bn((self.phi[ind2]-self.phi[ind1])/Vt)

    def dpJn(self,ind1,ind2,direction):
        if direction == '+': # derivative for ind1
            return self.n[ind1]*dBn((self.phi[ind1]-self.phi[ind2])/Vt)/Vt+self.n[ind2]*dBn((self.phi[ind2]-self.phi[ind1])/Vt)/Vt
        elif direction == '-': # derivative for ind2
            return -self.n[ind1]*dBn((self.phi[ind1]-self.phi[ind2])/Vt)/Vt-self.n[ind2]*dBn((self.phi[ind2]-self.phi[ind1])/Vt)/Vt
        else : print("wrong direction!")

    def dnJn(self,ind1,ind2,direction):
        if direction == '+':
            return Bn((self.phi[ind1]-self.phi[ind2])/Vt)
        elif direction == '-':
            return -Bn((self.phi[ind2]-self.phi[ind1])/Vt)
        else : print("wrong direction!")

    def update(self,phi,n):
        self.Rp = np.zeros(N)
        self.Rn = np.zeros(N)

        self.dpRp = np.zeros((N,N))
        self.dpRn = np.zeros((N,N))
        self.dnRp = np.zeros((N,N))
        self.dnRn = np.zeros((N,N))
   
        self.R = np.zeros(2*N)
        self.J = np.zeros((2*N,2*N))

        self.phi = copy.deepcopy(phi)
        self.n = copy.deepcopy(n)

    def make_RJ(self,Vg,VD):
                
	############
	## Vortex ##
	############

        ## left bottom cornor ##
        c = 0; u = Nx; r = 1

        self.Rp[c] = e1*(-self.phi[c]+0.5*self.phi[r]+0.5*self.phi[u])
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = 0.5*(self.Jn(r,c)+self.Jn(u,c))
        self.dpRp[c,c] = -e1; self.dpRp[c,r] = 0.5*e1; self.dpRp[c,u] = 0.5*e1
        #self.dpRn[c,c] = 0.5*(self.dpJn(r,c,'-')+self.dpJn(u,c,'-')); self.dpRn[c,r] = 0.5*self.dpJn(r,c,'+'); self.dpRn[c,u] = 0.5*self.dpJn(u,c,'+')
        self.dnRn[c,c] = 1. #0.5*(self.dnJn(r,c,'-')+self.dnJn(u,c,'-')); self.dnRn[c,r] = 0.5*self.dnJn(r,c,'+'); self.dnRn[c,u] = 0.5*self.dnJn(u,c,'+')

        ## right bottom cornor ##
        c = Nx-1; u = 2*Nx-1; l = Nx-2

        self.Rp[c] = e1*(-self.phi[c]+0.5*self.phi[l]+0.5*self.phi[u])
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = 0.5*(-self.Jn(c,l)+self.Jn(u,c))
        self.dpRp[c,c] = -e1; self.dpRp[c,l] = 0.5*e1; self.dpRp[c,u] = 0.5*e1
        #self.dpRn[c,c] = 0.5*(-self.dpJn(c,l,'+')+self.dpJn(u,c,'-')); self.dpRn[c,l] = -0.5*self.dpJn(c,l,'-'); self.dpRn[c,u] = 0.5*self.dpJn(u,c,'+')
        self.dnRn[c,c] = 1. #0.5*(-self.dnJn(c,l,'+')+self.dnJn(u,c,'-')); self.dnRn[c,l] = -0.5*self.dnJn(c,l,'-'); self.dnRn[c,u] = 0.5*self.dnJn(u,c,'+')

        ## left top cornor ##
        c = (Ny-1)*Nx; d = (Ny-2)*Nx; r = (Ny-1)*Nx+1

        self.Rp[c] = e1*(-self.phi[c]+0.5*self.phi[r]+0.5*self.phi[d])
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = 0.5*(self.Jn(r,c)-self.Jn(c,d))
        self.dpRp[c,c] = -e1; self.dpRp[c,r] = 0.5*e1; self.dpRp[c,d] = 0.5*e1
        #self.dpRn[c,c] = 0.5*(self.dpJn(r,c,'-')-self.dpJn(c,d,'+')); self.dpRn[c,r] = 0.5*self.dpJn(r,c,'+'); self.dpRn[c,d] = -0.5*self.dpJn(c,d,'-')
        self.dnRn[c,c] = 1. #0.5*(self.dnJn(r,c,'-')-self.dnJn(c,d,'+')); self.dnRn[c,r] = 0.5*self.dnJn(r,c,'+'); self.dnRn[c,d] = -0.5*self.dnJn(c,d,'-')

        ## right top cornor ##
        c = N-1; d = (Ny-1)*Nx-1; l = N-2

        self.Rp[c] = e1*(-self.phi[c]+0.5*self.phi[l]+0.5*self.phi[d])
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = 0.5*(-self.Jn(c,l)-self.Jn(c,d))
        self.dpRp[c,c] = -e1; self.dpRp[c,l]= 0.5*e1; self.dpRp[c,d] = 0.5*e1
        #self.dpRn[c,c] = 0.5*(-self.dpJn(c,l,'+')-self.dpJn(c,d,'+')); self.dpRn[c,l]= -0.5*self.dpJn(c,l,'-'); self.dpRn[c,d] = -0.5*self.dpJn(c,d,'-')
        self.dnRn[c,c] = 1. #0.5*(-self.dnJn(c,l,'+')-self.dnJn(c,d,'+')); self.dnRn[c,l]= -0.5*self.dnJn(c,l,'-'); self.dnRn[c,d] = -0.5*self.dnJn(c,d,'-')

	#######################
	## Edge at Interface ##
	#######################

        ## left bottom edge ##
        c = (My1-1)*Nx; u = My1*Nx; d = (My1-2)*Nx; r = (My1-1)*Nx+1

        self.Rp[c] = -(e1+e2)*self.phi[c]+(e1+e2)/2.*self.phi[r]+0.5*e1*self.phi[d]+0.5*e2*self.phi[u]
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = self.Jn(r,c)+0.5*self.Jn(u,c)-0.5*self.Jn(c,d)

        self.dpRp[c,c] = -(e1+e2); self.dpRp[c,r] = (e1+e2)/2.; self.dpRp[c,d] = 0.5*e1; self.dpRp[c,u] = 0.5*e2

        #self.dpRn[c,c] = self.dpJn(r,c,'-')+0.5*self.dpJn(u,c,'-')-0.5*self.dpJn(c,d,'+'); 
        #self.dpRn[c,r] = self.dpJn(r,c,'+'); self.dpRn[c,d] = -0.5*self.dpJn(c,d,'-'); self.dpRn[c,u] = 0.5*self.dpJn(u,c,'+') 

        self.dnRn[c,c] = 1. #self.dnJn(r,c,'-')+0.5*self.dnJn(u,c,'-')-0.5*self.dnJn(c,d,'+'); 
        #self.dnRn[c,r] = self.dnJn(r,c,'+'); self.dnRn[c,d] = -0.5*self.dnJn(c,d,'-'); self.dnRn[c,u] = 0.5*self.dnJn(u,c,'+') 
        
        ## right bottom edge ##
        c = My1*Nx-1; u = (My1+1)*Nx-1; d = (My1-1)*Nx-1; l = My1*Nx-2

        self.Rp[c] = -(e1+e2)*self.phi[c]+(e1+e2)/2.*self.phi[l]+0.5*e1*self.phi[d]+0.5*e2*self.phi[u]
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = -self.Jn(c,l)-0.5*self.Jn(c,d)+0.5*self.Jn(u,c)

        self.dpRp[c,c] = -(e1+e2); self.dpRp[c,l] = (e1+e2)/2.; self.dpRp[c,d] = 0.5*e1; self.dpRp[c,u] = 0.5*e2

        #self.dpRn[c,c] = -self.dpJn(c,l,'+')-0.5*self.dpJn(c,d,'+')+0.5*self.dpJn(u,c,'-'); 
        #self.dpRn[c,l] = -self.dpJn(c,l,'-'); self.dpRn[c,d] = -0.5*self.dpJn(c,d,'-'); self.dpRn[c,u] = 0.5*self.dpJn(u,c,'+') 

        self.dnRn[c,c] = 1. #-self.dnJn(c,l,'+')-0.5*self.dnJn(c,d,'+')+0.5*self.dnJn(u,c,'-'); 
        #self.dnRn[c,l] = -self.dnJn(c,l,'-'); self.dnRn[c,d] = -0.5*self.dnJn(c,d,'-'); self.dnRn[c,u] = 0.5*self.dnJn(u,c,'+') 

        ## left top edge ##
        c = My2*Nx; u = (My2+1)*Nx; d = (My2-1)*Nx; r = My2*Nx+1

        self.Rp[c] = -(e1+e2)*self.phi[c]+(e1+e2)/2.*self.phi[r]+0.5*e2*self.phi[d]+0.5*e1*self.phi[u]
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = self.Jn(r,c)-0.5*self.Jn(c,d)+0.5*self.Jn(u,c)

        self.dpRp[c,c] = -(e1+e2); self.dpRp[c,r] = (e1+e2)/2.; self.dpRp[c,d] = 0.5*e2; self.dpRp[c,u] = 0.5*e1

        #self.dpRn[c,c] = self.dpJn(r,c,'-')-0.5*self.dpJn(c,d,'+')+0.5*self.dpJn(u,c,'-'); 
        #self.dpRn[c,r] = self.dpJn(r,c,'+'); self.dpRn[c,d] = -0.5*self.dpJn(c,d,'-'); self.dpRn[c,u] = 0.5*self.dpJn(u,c,'+') 

        self.dnRn[c,c] = 1. #self.dnJn(r,c,'-')-0.5*self.dnJn(c,d,'+')+0.5*self.dnJn(u,c,'-'); 
        #self.dnRn[c,r] = self.dnJn(r,c,'+'); self.dnRn[c,d] = -0.5*self.dnJn(c,d,'-'); self.dnRn[c,u] = 0.5*self.dnJn(u,c,'+') 

        ## right top edge ##
        c = (My2+1)*Nx-1; u = (My2+2)*Nx-1; d = My2*Nx-1; l = (My2+1)*Nx-2

        self.Rp[c] = -(e1+e2)*self.phi[c]+(e1+e2)/2.*self.phi[l]+0.5*e2*self.phi[d]+0.5*e1*self.phi[u]
        self.Rn[c] = self.n[c] 
        #self.Rn[c] = -self.Jn(c,l)-0.5*self.Jn(c,d)+0.5*self.Jn(u,c)

        self.dpRp[c,c] = -(e1+e2); self.dpRp[c,l] = (e1+e2)/2.; self.dpRp[c,d] = 0.5*e2; self.dpRp[c,u] = 0.5*e1

        #self.dpRn[c,c] = -self.dpJn(c,l,'+')-0.5*self.dpJn(c,d,'+')+0.5*self.dpJn(u,c,'-'); 
        #self.dpRn[c,l] = -self.dpJn(c,l,'-'); self.dpRn[c,d] = -0.5*self.dpJn(c,d,'-'); self.dpRn[c,u] = 0.5*self.dpJn(u,c,'+') 

        self.dnRn[c,c] = 1. #-self.dnJn(c,l,'+')-0.5*self.dnJn(c,d,'+')+0.5*self.dnJn(u,c,'-'); 
        #self.dnRn[c,l] = -self.dnJn(c,l,'-'); self.dnRn[c,d] = -0.5*self.dnJn(c,d,'-'); self.dnRn[c,u] = 0.5*self.dnJn(u,c,'+') 

	########################
	## Bottom Oxide Layer ##
	########################

        for i in range(0,My1*Nx):
            ix = i%Nx
            iy = i/Nx
            if ix==0 and not(iy==0 or iy==My1-1) : 
                self.Rp[i] = e1*(-2.*self.phi[i]+self.phi[i+1]+0.5*self.phi[i+Nx]+0.5*self.phi[i-Nx])
                self.Rn[i] = self.n[i] 
                #self.Rn[i] = self.Jn(i+1,i)+0.5*self.Jn(i+Nx,i)-0.5*self.Jn(i,i-Nx)

                self.dpRp[i,i] = -2.*e1; self.dpRp[i,i+1] = e1; self.dpRp[i,i+Nx] = 0.5*e1; self.dpRp[i,i-Nx] = 0.5*e1

                #self.dpRn[i,i] = self.dpJn(i+1,i,'-')+0.5*self.dpJn(i+Nx,i,'-')-0.5*self.dpJn(i,i-Nx,'+')
                #self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i+Nx] = 0.5*self.dpJn(i+Nx,i,'+'); self.dpRn[i,i-Nx] = -0.5*self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #self.dnJn(i+1,i,'-')+0.5*self.dnJn(i+Nx,i,'-')-0.5*self.dnJn(i,i-Nx,'+')
                #self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i+Nx] = 0.5*self.dnJn(i+Nx,i,'+'); self.dnRn[i,i-Nx] = -0.5*self.dnJn(i,i-Nx,'-')

            elif ix==Nx-1 and not(iy==0 or iy==My1-1) : 
                self.Rp[i] = e1*(-2.*self.phi[i]+self.phi[i-1]+0.5*self.phi[i+Nx]+0.5*self.phi[i-Nx])
                self.Rn[i] = self.n[i] 
                #self.Rn[i] = -self.Jn(i,i-1)+0.5*self.Jn(i+Nx,i)-0.5*self.Jn(i,i-Nx)

                self.dpRp[i,i] = -2.*e1; self.dpRp[i,i-1] = e1; self.dpRp[i,i+Nx] = 0.5*e1; self.dpRp[i,i-Nx] = 0.5*e1

                #self.dpRn[i,i] = -self.dpJn(i,i-1,'+')+0.5*self.dpJn(i+Nx,i,'-')-0.5*self.dpJn(i,i-Nx,'+') 
                #self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i+Nx] = 0.5*self.dpJn(i+Nx,i,'+'); self.dpRn[i,i-Nx] = -0.5*self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #-self.dnJn(i,i-1,'+')+0.5*self.dnJn(i+Nx,i,'-')-0.5*self.dnJn(i,i-Nx,'+') 
                #self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i+Nx] = 0.5*self.dnJn(i+Nx,i,'+'); self.dnRn[i,i-Nx] = -0.5*self.dnJn(i,i-Nx,'-')

            elif iy==0 and not(ix==0 or ix==Nx-1) :
                if ix<Mx1 or ix>Mx2 :
                    self.Rp[i] = e1*(-2.*self.phi[i]+self.phi[i+Nx]+0.5*self.phi[i-1]+0.5*self.phi[i+1])
                    self.Rn[i] = self.n[i] 
                    #self.Rn[i] = self.Jn(i+Nx,i)+0.5*self.Jn(i+1,i)-0.5*self.Jn(i,i-1)

                    self.dpRp[i,i] = -2.*e1; self.dpRp[i,i+Nx] = e1; self.dpRp[i,i-1] = 0.5*e1; self.dpRp[i,i+1] = 0.5*e1

                    #self.dpRn[i,i] = self.dpJn(i+Nx,i,'-')+0.5*self.dpJn(i+1,i,'-')-0.5*self.dpJn(i,i-1,'+')
                    #self.dpRn[i,i+Nx] = self.dpJn(i+Nx,i,'+'); self.dpRn[i,i+1] = 0.5*self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -0.5*self.dpJn(i,i-1,'-')

                    self.dnRn[i,i] = 1. #self.dnJn(i+Nx,i,'-')+0.5*self.dnJn(i+1,i,'-')-0.5*self.dnJn(i,i-1,'+')
                    #self.dnRn[i,i+Nx] = self.dnJn(i+Nx,i,'+'); self.dnRn[i,i+1] = 0.5*self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -0.5*self.dnJn(i,i-1,'-')
                else :
                    self.dpRp[i,i] = 1.; self.Rp[i] = self.phi[i] - (phi0+Vg)
                    self.dnRn[i,i] = 1.; self.Rn[i] = self.n[i]

            elif iy==My1-1 and not(ix==0 or ix==Nx-1) : 
                self.Rp[i] = -2.*(e1+e2)*self.phi[i]+(e1+e2)/2.*self.phi[i-1]+(e1+e2)/2.*self.phi[i+1]+e1*self.phi[i-Nx]+e2*self.phi[i+Nx]
                self.Rn[i] = self.n[i]
                #self.Rn[i] = self.Jn(i+Nx,i)+self.Jn(i+1,i)-self.Jn(i,i-Nx)-self.Jn(i,i-1)

                self.dpRp[i,i] = -2.*(e1+e2); self.dpRp[i,i-Nx] = e1; self.dpRp[i,i+Nx] = e2; self.dpRp[i,i-1] = 0.5*(e1+e2); self.dpRp[i,i+1] = 0.5*(e1+e2)

                #self.dpRn[i,i] = self.dpJn(i+Nx,i,'-')+self.dpJn(i+1,i,'-')-self.dpJn(i,i-1,'+')-self.dpJn(i,i-Nx,'+')
                #self.dpRn[i,i+Nx] = self.dpJn(i+Nx,i,'+'); self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i-Nx] = -self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #self.dnJn(i+Nx,i,'-')+self.dnJn(i+1,i,'-')-self.dnJn(i,i-1,'+')-self.dnJn(i,i-Nx,'+')
                #self.dnRn[i,i+Nx] = self.dnJn(i+Nx,i,'+'); self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i-Nx] = -self.dnJn(i,i-Nx,'-')

            elif not(ix==0 or ix==Nx-1) and not(iy==0 or iy==My1-1) :
                self.Rp[i] = e1*(-4.*self.phi[i]+self.phi[i-1]+self.phi[i+1]+self.phi[i-Nx]+self.phi[i+Nx])
                self.Rn[i] = self.n[i]
                #self.Rn[i] = self.Jn(i+Nx,i)+self.Jn(i+1,i)-self.Jn(i,i-Nx)-self.Jn(i,i-1)

                self.dpRp[i,i] = -4.*e1; self.dpRp[i,i-Nx] = e1; self.dpRp[i,i+Nx] = e1; self.dpRp[i,i-1] = e1; self.dpRp[i,i+1] = e1

                #self.dpRn[i,i] = self.dpJn(i+Nx,i,'-')+self.dpJn(i+1,i,'-')-self.dpJn(i,i-1,'+')-self.dpJn(i,i-Nx,'+')
                #self.dpRn[i,i+Nx] = self.dpJn(i+Nx,i,'+'); self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i-Nx] = -self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #self.dnJn(i+Nx,i,'-')+self.dnJn(i+1,i,'-')-self.dnJn(i,i-1,'+')-self.dnJn(i,i-Nx,'+')
                #self.dnRn[i,i+Nx] = self.dnJn(i+Nx,i,'+'); self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i-Nx] = -self.dnJn(i,i-Nx,'-')
        
	###################
	## Silicon Layer ##
	###################

        for i in range(My1*Nx,My2*Nx):
            ix = i%Nx
            iy = i/Nx
            if ix==0 : 
                self.dpRp[i,i] = 1.; self.Rp[i] = self.phi[i] - phi1
                self.dnRn[i,i] = 1.; self.Rn[i] = self.n[i] - np.abs(Np)
            elif ix==Nx-1 : 
                self.dpRp[i,i] = 1.; self.Rp[i] = self.phi[i] - (phi1+VD)
                self.dnRn[i,i] = 1.; self.Rn[i] = self.n[i] - np.abs(Np)
            else :
                if ix<Mx1 or ix>Mx2:
                    self.Rp[i] = e2*(-4.*self.phi[i]+self.phi[i-1]+self.phi[i+1]+self.phi[i-Nx]+self.phi[i+Nx]) - h*h*const.e*(Np+self.n[i])/const.epsilon_0
                else :
                    self.Rp[i] = e2*(-4.*self.phi[i]+self.phi[i-1]+self.phi[i+1]+self.phi[i-Nx]+self.phi[i+Nx]) - h*h*const.e*self.n[i]/const.epsilon_0

                self.Rn[i] = self.Jn(i+Nx,i)+self.Jn(i+1,i)-self.Jn(i,i-Nx)-self.Jn(i,i-1)

                self.dpRp[i,i] = -4.*e2; self.dpRp[i,i-Nx] = e2; self.dpRp[i,i+Nx] = e2; self.dpRp[i,i-1] = e2; self.dpRp[i,i+1] = e2
                self.dnRp[i,i] = -h*h*const.e/const.epsilon_0;

                self.dpRn[i,i] = self.dpJn(i+Nx,i,'-')+self.dpJn(i+1,i,'-')-self.dpJn(i,i-1,'+')-self.dpJn(i,i-Nx,'+')
                self.dpRn[i,i+Nx] = self.dpJn(i+Nx,i,'+'); self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i-Nx] = -self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = self.dnJn(i+Nx,i,'-')+self.dnJn(i+1,i,'-')-self.dnJn(i,i-1,'+')-self.dnJn(i,i-Nx,'+')
                self.dnRn[i,i+Nx] = self.dnJn(i+Nx,i,'+'); self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i-Nx] = -self.dnJn(i,i-Nx,'-')

	#####################
	## Top Oxide Layer ##
	#####################

        for i in range(My2*Nx,N):
            ix = i%Nx
            iy = i/Nx
            if ix==0 and not(iy==My2 or iy==Ny-1) : 
                self.Rp[i] = e1*(-2.*self.phi[i]+self.phi[i+1]+0.5*self.phi[i+Nx]+0.5*self.phi[i-Nx])
                self.Rn[i] = self.n[i] 
                #self.Rn[i] = self.Jn(i+1,i)+0.5*self.Jn(i+Nx,i)-0.5*self.Jn(i,i-Nx)

                self.dpRp[i,i] = -2.*e1; self.dpRp[i,i+1] = e1; self.dpRp[i,i+Nx] = 0.5*e1; self.dpRp[i,i-Nx] = 0.5*e1

                #self.dpRn[i,i] = self.dpJn(i+1,i,'-')+0.5*self.dpJn(i+Nx,i,'-')-0.5*self.dpJn(i,i-Nx,'+')
                #self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i+Nx] = 0.5*self.dpJn(i+Nx,i,'+'); self.dpRn[i,i-Nx] = -0.5*self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #self.dnJn(i+1,i,'-')+0.5*self.dnJn(i+Nx,i,'-')-0.5*self.dnJn(i,i-Nx,'+')
                #self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i+Nx] = 0.5*self.dnJn(i+Nx,i,'+'); self.dnRn[i,i-Nx] = -0.5*self.dnJn(i,i-Nx,'-')

            elif ix==Nx-1 and not(iy==My2 or iy==Ny-1) : 
                self.Rp[i] = e1*(-2.*self.phi[i]+self.phi[i-1]+0.5*self.phi[i+Nx]+0.5*self.phi[i-Nx])
                self.Rn[i] = self.n[i] 
                #self.Rn[i] = -self.Jn(i,i-1)+0.5*self.Jn(i+Nx,i)-0.5*self.Jn(i,i-Nx)

                self.dpRp[i,i] = -2.*e1; self.dpRp[i,i-1] = e1; self.dpRp[i,i+Nx] = 0.5*e1; self.dpRp[i,i-Nx] = 0.5*e1

                #self.dpRn[i,i] = -self.dpJn(i,i-1,'+')+0.5*self.dpJn(i+Nx,i,'-')-0.5*self.dpJn(i,i-Nx,'+') 
                #self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i+Nx] = 0.5*self.dpJn(i+Nx,i,'+'); self.dpRn[i,i-Nx] = -0.5*self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #-self.dnJn(i,i-1,'+')+0.5*self.dnJn(i+Nx,i,'-')-0.5*self.dnJn(i,i-Nx,'+') 
                #self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i+Nx] = 0.5*self.dnJn(i+Nx,i,'+'); self.dnRn[i,i-Nx] = -0.5*self.dnJn(i,i-Nx,'-')

            elif iy==My2 and not(ix==0 or ix==Nx-1) :
                self.Rp[i] = -2.*(e1+e2)*self.phi[i]+(e1+e2)/2.*self.phi[i-1]+(e1+e2)/2.*self.phi[i+1]+e2*self.phi[i-Nx]+e1*self.phi[i+Nx]
		self.Rn[i] = self.n[i]
                #self.Rn[i] = self.Jn(i+Nx,i)+self.Jn(i+1,i)-self.Jn(i,i-Nx)-self.Jn(i,i-1)

                self.dpRp[i,i] = -2.*(e1+e2); self.dpRp[i,i-Nx] = e2; self.dpRp[i,i+Nx] = e1; self.dpRp[i,i-1] = 0.5*(e1+e2); self.dpRp[i,i+1] = 0.5*(e1+e2)

                #self.dpRn[i,i] = self.dpJn(i+Nx,i,'-')+self.dpJn(i+1,i,'-')-self.dpJn(i,i-1,'+')-self.dpJn(i,i-Nx,'+')
                #self.dpRn[i,i+Nx] = self.dpJn(i+Nx,i,'+'); self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i-Nx] = -self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #self.dnJn(i+Nx,i,'-')+self.dnJn(i+1,i,'-')-self.dnJn(i,i-1,'+')-self.dnJn(i,i-Nx,'+')
                #self.dnRn[i,i+Nx] = self.dnJn(i+Nx,i,'+'); self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i-Nx] = -self.dnJn(i,i-Nx,'-')

            elif iy==Ny-1 and not(ix==0 or ix==Nx-1) :
                if ix<Mx1 or ix>Mx2 :
                    self.Rp[i] = e1*(-2.*self.phi[i]+self.phi[i-Nx]+0.5*self.phi[i-1]+0.5*self.phi[i+1])
                    self.Rn[i] = self.n[i] 
                    #self.Rn[i] = -self.Jn(i,i-Nx)+0.5*self.Jn(i+1,i)-0.5*self.Jn(i,i-1)

                    self.dpRp[i,i] = -2.*e1; self.dpRp[i,i-Nx] = e1; self.dpRp[i,i-1] = 0.5*e1; self.dpRp[i,i+1] = 0.5*e1

                    #self.dpRn[i,i] = -self.dpJn(i,i-Nx,'+')+0.5*self.dpJn(i+1,i,'-')-0.5*self.dpJn(i,i-1,'+')
                    #self.dpRn[i,i-Nx] = -self.dpJn(i,i-Nx,'-'); self.dpRn[i,i+1] = 0.5*self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -0.5*self.dpJn(i,i-1,'-')

                    self.dnRn[i,i] = 1. #-self.dnJn(i,i-Nx,'+')+0.5*self.dnJn(i+1,i,'-')-0.5*self.dnJn(i,i-1,'+')
                    #self.dnRn[i,i-Nx] = -self.dnJn(i,i-Nx,'-'); self.dnRn[i,i+1] = 0.5*self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -0.5*self.dnJn(i,i-1,'-')
                else :
                    self.dpRp[i,i] = 1.; self.Rp[i] = self.phi[i] - (phi0+Vg)
                    self.dnRn[i,i] = 1.; self.Rn[i] = self.n[i]

            elif not(ix==0 or ix==Nx-1) and not(iy==My2 or iy==Ny-1) :
                self.Rp[i] = e1*(-4.*self.phi[i]+self.phi[i-1]+self.phi[i+1]+self.phi[i-Nx]+self.phi[i+Nx])
                self.Rn[i] = self.n[i]
                #self.Rn[i] = self.Jn(i+Nx,i)+self.Jn(i+1,i)-self.Jn(i,i-Nx)-self.Jn(i,i-1)

                self.dpRp[i,i] = -4.*e1; self.dpRp[i,i-Nx] = e1; self.dpRp[i,i+Nx] = e1; self.dpRp[i,i-1] = e1; self.dpRp[i,i+1] = e1

                #self.dpRn[i,i] = self.dpJn(i+Nx,i,'-')+self.dpJn(i+1,i,'-')-self.dpJn(i,i-1,'+')-self.dpJn(i,i-Nx,'+')
                #self.dpRn[i,i+Nx] = self.dpJn(i+Nx,i,'+'); self.dpRn[i,i+1] = self.dpJn(i+1,i,'+'); self.dpRn[i,i-1] = -self.dpJn(i,i-1,'-'); self.dpRn[i,i-Nx] = -self.dpJn(i,i-Nx,'-')

                self.dnRn[i,i] = 1. #self.dnJn(i+Nx,i,'-')+self.dnJn(i+1,i,'-')-self.dnJn(i,i-1,'+')-self.dnJn(i,i-Nx,'+')
                #self.dnRn[i,i+Nx] = self.dnJn(i+Nx,i,'+'); self.dnRn[i,i+1] = self.dnJn(i+1,i,'+'); self.dnRn[i,i-1] = -self.dnJn(i,i-1,'-'); self.dnRn[i,i-Nx] = -self.dnJn(i,i-Nx,'-')

        self.J[0::2,0::2] = self.dpRp
        self.J[0::2,1::2] = self.dnRp
        self.J[1::2,0::2] = self.dpRn
        self.J[1::2,1::2] = self.dnRn

        self.R[0::2] = self.Rp
        self.R[1::2] = self.Rn

        self.J = np.matmul(self.J,self.Cmat)
        Rvec = 1./np.sum(np.abs(self.J),axis=1); Rmat = np.diag(Rvec)
        self.J = np.matmul(Rmat,self.J)
        self.R = np.matmul(Rmat,self.R)
   
    def make_Current(self):
        Crnt1 = np.zeros((Ny,Nx))
        Crnt2 = np.zeros((Ny,Nx))

	for i in range(N):
	    ix = i%Nx
	    iy = i/Nx
	    if not(ix==0 or ix==Nx-1) and not(iy==0 or iy==Ny-1):
	        Crnt1[iy,ix] = const.e*Dn*self.Jn(i,i-1)*h
		#Crnt2[iy,ix] = const.e*Dn*self.Jn(i,i-Ny)*h

        return Crnt1,Crnt2
 
def Newton_DD(Drift_Diffusion,phin,dx_list,Vg,VD):
    dx = 1.;
    while dx > 1e-6:
        Drift_Diffusion.make_RJ(Vg,VD)
        R = Drift_Diffusion.R
        J = Drift_Diffusion.J

        #plt.imshow(J); plt.show()
        #plt.imshow(R[0::2].reshape((Ny,Nx))); plt.show()
        #plt.imshow(R[1::2].reshape((Ny,Nx))); plt.show()

        dphin = spsolve(csr_matrix(J),-R)
        phin += np.matmul(Drift_Diffusion.Cmat,dphin)
        Drift_Diffusion.update(phin[0::2],phin[1::2])

        dx = np.linalg.norm(dphin); print(dx)
        dx_list.append(dx)

    return phin

def Current(Drift_Diffusion):
    I = 0
    for yi in range(My1,My2):
        I += const.e*Dn*Drift_Diffusion.Jn(yi*Nx+Nx-1,yi*Nx+Nx-2); # Drain Current
        #I += const.e*Dn*Drift_Diffusion.Jn(yi*Nx+1,yi*Nx); # Source Current

    return I*h

Vg = 1.1; VD = 1.1; dx_list = []

#Prob = Poisson_eq_2D(Vg)
#phi = spsolve(csr_matrix(Prob.H.csr_form(),shape=(N,N)),Prob.b) ### initial vector 1
#phi = Newton(Prob,phi,dx_list)
#n = Elec_density(phi)

phi = np.zeros(N)
n = np.ones(N)

phin = np.zeros(2*N)
phin[0::2] = phi; phin[1::2] = n
Prob = Drift_Diffusion_2D(phin[0::2],phin[1::2])
'''
phin = Newton_DD(Prob,phin,dx_list,Vg,VD)
phi = phin[0::2]; n = phin[1::2]

fig = plt.figure(figsize=(15,5))
ax1 = fig.add_subplot(121,projection='3d')
ax2 = fig.add_subplot(122,projection='3d')

phi = phi.reshape((Ny,Nx))
n = n.reshape((Ny,Nx))

X,Y = np.meshgrid(np.arange(0,Nx)*h,np.arange(0,Ny)*h)
surf1 = ax1.plot_surface(X,Y,phi,cmap=cm.coolwarm,ccount=150)
surf2 = ax2.plot_surface(X,Y,n,edgecolors='r',alpha=0)

ax1.set_xlabel('x',fontsize=18)
ax1.set_ylabel('y',fontsize=18)
ax1.set_zlabel('$\phi(x,y)$',fontsize=18)
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)
ax1.text(0,0,0.45,'$V_G=$'+str(Vg),fontsize=18)

ax2.set_xlabel('x',fontsize=18)
ax2.set_ylabel('y',fontsize=18)
ax2.set_zlabel('$n(x,y)$',fontsize=18)
ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0)

fig.colorbar(surf1,ax=ax1)
fig.tight_layout()
plt.show()

Crnt1,Crnt2 = Prob.make_Current()

fig = plt.figure(figsize=(10,10))
plt.quiver(X,Y,Crnt1,Crnt2,angles='xy',pivot='middle',linewidths=np.linspace(0,1.,X.size))

plt.show()

'''
fig = plt.figure(figsize=(10,10))

VD_list = np.linspace(0,1.1,11)
Vg_list = np.linspace(0,1.1,3)
I_list = []
for Vg in Vg_list:
    for VD in VD_list:
        phin = Newton_DD(Prob,phin,dx_list,Vg,VD)
        I = Current(Prob)
        I_list.append(I)

    plt.plot(VD_list,I_list,'-',c=plt.cm.coolwarm(Vg/1.1),label='$V_G=$'+str(format(Vg,'1.1f'))+'V')
    I_list = []

plt.xlabel('$V_D$',fontsize=18)
plt.ylabel('$I_D$',fontsize=18)
plt.legend(loc='best',frameon=False,fontsize=15)

plt.show()

