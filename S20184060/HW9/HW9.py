import numpy as np
from scipy.linalg import solve
import scipy.constants as const 
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

Nx = 90
Ny = 50
N = Nx*Ny

class Poisson_eq_2D:
    def __init__(self):
        self.H = np.zeros((N,N))
        self.b = np.zeros(N)
        
        self.H[0,0] = -1.; self.H[0,1] = 0.5; self.H[0,Nx] = 0.5
        self.H[Nx-1,Nx-1] = -1.; self.H[Nx-1,Nx-2] = 0.5; self.H[Nx-1,2*Nx-1] = 0.5
        self.H[(Ny-1)*Nx,(Ny-1)*Nx] = -1.; self.H[(Ny-1)*Nx,(Ny-1)*Nx+1] = 0.5; self.H[(Ny-1)*Nx,(Ny-2)*Nx] = 0.5
        self.H[N-1,N-1] =  -1.; self.H[N-1,N-2] = 0.5; self.H[N-1,(Ny-1)*Nx-1] = 0.5

        for i in range(N):
            ix = i%Nx
            iy = i/Nx
            if ix==0 and not(iy==0 or iy==Ny-1) : 
                self.H[i,i] = -2.; self.H[i,i+1] = 1.; self.H[i,i+Nx] = 0.5; self.H[i,i-Nx] = 0.5
            elif ix==Nx-1 and not(iy==0 or iy==Ny-1) : 
                self.H[i,i] = -2.; self.H[i,i-1] = 1.; self.H[i,i+Nx] = 0.5; self.H[i,i-Nx] = 0.5
            elif iy==0 and not(ix==0 or ix==Nx-1) :
                self.H[i,i] = -2.; self.H[i,i+Nx] = 1.; self.H[i,i-1] = 0.5; self.H[i,i+1] = 0.5
            elif iy==Ny-1 and not(ix==0 or ix==Nx-1) :
                self.H[i,i] = -2.; self.H[i,i-Nx] = 1.; self.H[i,i-1] = 0.5; self.H[i,i+1] = 0.5
            elif not(ix==0 or ix==Nx-1) and not(iy==0 or iy==Ny-1) :
                self.H[i,i] = -4.; self.H[i,i-1] = 1.; self.H[i,i+1] = 1.; self.H[i,i-Nx] = 1.; self.H[i,i+Nx] = 1.

    def make_prob(self,BC_list):
        for ix,iy,bi in BC_list:
            i = iy*Nx+ix
            H_row = np.zeros(N); H_row[i] = 1.
            self.H[i,:] = H_row
            self.b[i] = bi

def make_BC(BC_type):
    BC_list = []

    if BC_type==0:
        for ix in range(Nx):
            BC_list.append((ix,0,0.))
        for ix in range((Nx*2)/9):
            BC_list.append((ix,Ny-1,0.))
            BC_list.append((Nx-1-ix,Ny-1,1.))
    elif BC_type==1:
        for ix in range(Nx):
            BC_list.append((ix,0,0.))
        for ix in range((Nx*2)/9):
            BC_list.append((ix,Ny-1,1.))
            BC_list.append((Nx-1-ix,Ny-1,0.))
    elif BC_type==2:
        for ix in range(Nx):
            BC_list.append((ix,0,1.))
        for ix in range((Nx*2)/9):
            BC_list.append((ix,Ny-1,0.))
            BC_list.append((Nx-1-ix,Ny-1,0.))
    elif BC_type==3:
        for ix in range(Nx):
            BC_list.append((ix,0,1.))
        for ix in range((Nx*2)/9):
            BC_list.append((ix,Ny-1,1.))
            BC_list.append((Nx-1-ix,Ny-1,1.))

    return BC_list

Prob = Poisson_eq_2D()
BC_list = make_BC(0)
Prob.make_prob(BC_list)

phi = solve(Prob.H,Prob.b)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111,projection='3d')

phi = phi.reshape((Ny,Nx))
X,Y = np.meshgrid(np.arange(0,Nx),np.arange(0,Ny))
surf = ax.plot_surface(X,Y,phi,cmap=cm.coolwarm)

ax.set_xlabel('x',fontsize=18)
ax.set_ylabel('y',fontsize=18)
ax.set_zlabel('$\phi(x,y)$',fontsize=18)
ax.set_xlim(0,Nx-1)
ax.set_ylim(0,Ny-1)
ax.set_zlim(0,1)

fig.colorbar(surf)
plt.show()
