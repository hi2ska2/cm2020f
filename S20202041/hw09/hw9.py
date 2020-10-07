import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

nx=9
ny=5

red = 1
blue = 0

problem = 4
phi = np.zeros(nx*ny)
if(problem==1):
	phi[nx*0+0] = blue
	phi[nx*0+1] = blue
	phi[nx*0+nx-2] = red
	phi[nx*0+nx-1] = red
	phi[nx*(ny-1):nx*ny] = blue
elif(problem==2) :
	phi[nx*0+0] = red
	phi[nx*0+1] = red
	phi[nx*0+nx-2] = blue
	phi[nx*0+nx-1] = blue
	phi[nx*(ny-1):nx*ny] = blue
elif(problem==3):
	phi[nx*0+0] = blue
	phi[nx*0+1] = blue
	phi[nx*0+nx-2] = blue
	phi[nx*0+nx-1] = blue
	phi[nx*(ny-1):nx*ny] = red
else:
	phi[nx*0+0] = red
	phi[nx*0+1] = red
	phi[nx*0+nx-2] = red
	phi[nx*0+nx-1] = red
	phi[nx*(ny-1):nx*ny] = red



H = np.zeros((nx*ny,nx*ny))
#blue and red
np.fill_diagonal(H[0:nx*ny],1)
np.fill_diagonal(H[0:nx*(ny-1)],0)
H[0,0] = 1
H[1,1] = 1
H[nx-2,nx-2] = 1
H[nx-1,nx-1] = 1


#black and white
#upper boundary black
for i in range(5):
	H[2+i,2+i] = -2; H[2+i,2+i+1]=0.5 ; H[2+i,2+i-1] = 0.5; H[2+i,2+i+nx] = 1;
#left boundary black
for i in range(3):
	H[nx+i*nx,nx+i*nx] = -2 ; H[nx+i*nx,nx+nx+i*nx] =0.5; H[nx+i*nx,i*nx]=0.5; H[nx+i*nx,nx+1+i*nx]=1;
#right boundary black
for i in range(3):
	H[2*nx-1+i*nx,2*nx-1+i*nx] = -2 ; H[2*nx-1+i*nx,2*nx-1+nx+i*nx]=0.5;H[2*nx-1+i*nx,2*nx-1-nx+i*nx]=0.5; H[2*nx-1+i*nx,2*nx-1-1+i*nx] = 1;
#all whites
for i in range(7):
	for j in range(3):
		H[nx+1+i+nx*j,nx+1+i+nx*j] = -4; H[nx+1+i+nx*j,nx+1+1+i+nx*j] = 1;H[nx+1+i+nx*j,nx+1-1+i+nx*j] = 1; H[nx+1+i+nx*j,nx+1+nx+i+nx*j] = 1;H[nx+1+i+nx*j,nx+1-nx+i+nx*j] = 1;


real_phi = slin.solve(H,phi)
Z = np.reshape(real_phi,(-1,nx))
x = range(nx)
y = range(ny)
X,Y = np.meshgrid(x,y)

Z = np.log(Z)


norm = plt.Normalize(Z.min(), Z.max())
colors = cm.viridis(norm(Z))
rcount, ccount, _ = colors.shape

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.zaxis.set_rotate_label(False)
ax.set_xlabel("x",fontsize=12)
ax.set_ylabel("y",fontsize=12)
ax.set_zlabel(r"$\log \phi(x,y)$",fontsize=12,rotation=90)
#ax.set_zscale('log')
ax.tick_params(labelsize=12)
surf = ax.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount,
                       facecolors=colors, shade=False)
surf.set_facecolor((0,0,0,0))
plt.show()


'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.zaxis.set_rotate_label(False)
ax.set_xlabel("x",fontsize=12)
ax.set_ylabel("y",fontsize=12)
ax.set_zlabel(r"$\log \phi(x,y)$",fontsize=12,rotation=90)
#ax.set_zscale('log')
ax.tick_params(labelsize=12)
ax.plot_wireframe(X, Y, Z)
plt.show()
'''

'''
plt.plot(x, e_density,label=r'$V_{gate}=0.3V$',c='r',marker='o',lw=0,ms=10,mfc='none')
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$n(r)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
'''
