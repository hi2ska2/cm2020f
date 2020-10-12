import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm


nx = 121
ny = 29

#coefficients
Nc = 2.86*1e+25
kT = sc.k*300
ni = Nc*np.exp(-0.56*sc.e/kT)
Nacc = 1e+26
q = sc.e 
e1 = 3.9*sc.epsilon_0
e2 = 11.7*sc.epsilon_0

dx = 30*1e-9/(nx-1)


GV = 1.1

phi = np.full(nx*ny,0.33374)
dphi_bc = np.zeros((nx*ny,nx*ny))

#setting phi_bc and corresponding H
H = np.zeros((nx*ny,nx*ny))

interx1 = int(nx/3)
interx2 = int(2*nx/3)
intery1 = int(ny/7)
intery2 = int(6*ny/7)

phi_bc = np.zeros(nx*ny)
phi_bc[interx1:interx2+1] = 0.33374 + GV
phi_bc[nx*(ny-1)+interx1:nx*(ny-1)+interx2+1] = 0.33374+GV

np.fill_diagonal(H[0:nx*(ny-1)+interx2+1],1)
np.fill_diagonal(H[0:nx*(ny-1)+interx1],0)
np.fill_diagonal(H[0:interx2+1],1)
np.fill_diagonal(H[0:interx1],0)

for i in range(intery2-intery1-1) :
	phi_bc[nx*int(intery1+1 + i)] = np.arcsinh(Nacc/(2*ni))*kT/q
	phi_bc[nx*int(intery1+1 + i + 1) -1] = np.arcsinh(Nacc/(2*ni))*kT/q
	H[nx*int(intery1+1+i),nx*int(intery1+1+i)] = 1
	H[nx*int(intery1+1+i+1)-1,nx*int(intery1+1+i+1)-1] = 1


#Neumann Boundary Conditions
#All corners
H[0,0] = -e1; H[0,1] = 0.5*e1 ; H[0,nx] = 0.5*e1;
H[nx-1,nx-1] = -e1; H[nx-1,nx-2] =0.5*e1; H[nx-1,2*nx-1] = 0.5*e1;
H[nx*(ny-1),nx*(ny-1)] = -e1; H[nx*(ny-1),nx*(ny-1)+1] = 0.5*e1; H[nx*(ny-1),nx*(ny-2)] = 0.5*e1;
H[nx*ny-1,nx*ny-1] = -e1; H[nx*ny-1,nx*ny-2] = 0.5*e1; H[nx*ny-1,nx*(ny-1)-1] = 0.5*e1;

#Upper&below boundary
for i in range(interx1-1):
	H[1+i,1+i] = -2*e1; H[1+i,1+i-1] = 0.5*e1; H[1+i,1+i+1] = 0.5*e1; H[1+i,1+i+nx] = e1;
	ru = interx2 + i +1
	H[ru,ru] = -2*e1; H[ru,ru-1] = 0.5*e1; H[ru,ru+1] = 0.5*e1; H[ru,ru+nx] = e1;
	leftb = nx*(ny-1)+1+i
	H[leftb,leftb] = -2*e1; H[leftb,leftb+1] = 0.5*e1; H[leftb,leftb-1]=0.5*e1; H[leftb,leftb-nx] = e1;
	rightb = nx*(ny-1)+1+i+interx2
	H[rightb,rightb] = -2*e1; H[rightb,rightb+1] = 0.5*e1; H[rightb,rightb-1]=0.5*e1; H[rightb,rightb-nx] = e1;

#left&right boundary
for i in range(intery1-1):
	lu = (i+1)*nx
	H[lu,lu] = -2*e1; H[lu,lu-nx] = 0.5*e1; H[lu,lu+nx] = 0.5*e1; H[lu,lu+1] = e1;
	rightu = (i+2)*nx-1
	H[rightu,rightu] = -2*e1; H[rightu,rightu-nx] = 0.5*e1; H[rightu,rightu+nx] = 0.5*e1; H[rightu,rightu-1] = e1;
	ld = nx*(i+intery2+1)
	rd = nx*(i+intery2+2)-1
	H[ld,ld] = -2*e1;H[ld,ld-nx] = 0.5*e1;H[ld,ld+nx]=0.5*e1; H[ld,ld+1]=e1;
	H[rd,rd] = -2*e1;H[rd,rd-nx] = 0.5*e1;H[rd,rd+nx]=0.5*e1; H[rd,rd+1]=e1;
	
#interfaces
lu = nx*intery1; ru = nx*intery1+nx-1;
ld = nx*intery2; rd = nx*intery2+nx-1;

H[lu,lu] = -(e1+e2); H[lu,lu-nx] = 0.5*e1; H[lu,lu+nx] = 0.5*e2; H[lu,lu+1] = (e1+e2)/2;
H[ru,ru] = -(e1+e2); H[ru,ru-nx] = 0.5*e1; H[ru,ru+nx] = 0.5*e2; H[ru,ru-1] = (e1+e2)/2;
H[ld,ld]= -(e1+e2);H[ld,ld-nx] = 0.5*e2;H[ld,ld+nx]=0.5*e1; H[ld,ld+1]=(e1+e2)/2;
H[rd,rd]= -(e1+e2);H[rd,rd-nx] = 0.5*e2;H[rd,rd+nx]=0.5*e1; H[rd,rd+1]=(e1+e2)/2;

#middle part
for j in range(ny-2):
	for i in range(nx-2):
		centeru = nx*(j+1)+i+1
		if(j+1<intery1) : H[centeru,centeru] = -4*e1; H[centeru,centeru-1] = e1; H[centeru,centeru+1] = e1; H[centeru,centeru+nx]=e1;H[centeru,centeru-nx] = e1;
		elif(j+1==intery1) : H[centeru,centeru] = -2*(e1+e2); H[centeru,centeru-1] = (e1+e2)/2; H[centeru,centeru+1] = (e1+e2)/2; H[centeru,centeru+nx]=e2;H[centeru,centeru-nx] = e1;
		elif((j+1>intery1)&(j+1<intery2)) : 
			H[centeru,centeru] = -4*e2; H[centeru,centeru-1] = e2; H[centeru,centeru+1] = e2; H[centeru,centeru+nx]=e2; H[centeru,centeru-nx] = e2;
			if((i+1>=interx1)&(i+1<=interx2)):
				phi_bc[centeru] = dx*dx*q*2*ni*np.sinh(q*phi[centeru]/kT)
				dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*q/kT
			else : 
				phi_bc[centeru] = dx*dx*q*(Nacc+2*ni*np.sinh(q*phi[centeru]/kT))
				dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*q/kT
		elif(j+1==intery2) : H[centeru,centeru] = -2*(e1+e2); H[centeru,centeru-1] = (e1+e2)/2; H[centeru,centeru+1] = (e1+e2)/2; H[centeru,centeru+nx]=e1;H[centeru,centeru-nx] = e2;
		else : H[centeru,centeru] = -4*e1; H[centeru,centeru-1] = e1; H[centeru,centeru+1] = e1; H[centeru,centeru+nx]=e1;H[centeru,centeru-nx] = e1;
			
print("for all done")

for i in range(nx*ny):
	tmp = 0
	for j in range(ny*nx):
		if(H[i,j]!=0) : tmp=tmp+1
	if(tmp==0): 
		yline = int(i/nx)
		xline = int(i%nx)
		print("(",yline,xline,")")


tmp = 200
min_phi = 100 
itr = 0

centerphi = np.zeros(39)

while abs(tmp-min_phi)>1e-5 :
	print("itr : ",itr+1,"abs(diff) : ",tmp-min_phi)
	min_phi = np.min(phi)
	residue = np.matmul(H,phi) - phi_bc
	Jaco = H -dphi_bc
	phi = phi - slin.solve(Jaco,residue)
	tmp = np.min(phi)
	centerphi[itr] = abs(tmp-min_phi)
	for j in range(ny-2):
		for i in range(nx-2):
			centeru = nx*(j+1)+i+1
			if((j+1>intery1)&(j+1<intery2)) : 
				if((i+1>=interx1)&(i+1<=interx2)):
					phi_bc[centeru] = dx*dx*q*2*ni*np.sinh(q*phi[centeru]/kT)
					dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*q/kT
				else : 
					phi_bc[centeru] = dx*dx*q*(Nacc+2*ni*np.sinh(q*phi[centeru]/kT))
					dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*q/kT
	itr = itr+1			

plt.plot(range(itr), centerphi[1:],label=r'$V_{gate}=1.1$',c='r',marker='o',lw=0,ms=10,mfc='none')
plt.xlabel("Iterations",fontsize = 18)
plt.ylabel(r'$|\phi^{(n+1)}_{min}-\phi^{(n)}_{min}|$',fontsize = 18)
plt.legend(fontsize=18)
plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()

	
Z = np.reshape(phi,(-1,nx))
x = np.linspace(0,30e-9,nx,endpoint=True)
y = np.linspace(0,7e-9,ny,endpoint=True)
X,Y = np.meshgrid(x,y)

#Z = np.log(Z)

norm = plt.Normalize(Z.min(), Z.max())
colors = cm.viridis(norm(Z))
rcount, ccount, _ = colors.shape

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.zaxis.set_rotate_label(False)
ax.set_xlabel("x",fontsize=12)
ax.set_ylabel("y",fontsize=12)
ax.set_zlabel(r"$\phi(x,y)$",fontsize=12,rotation=90)
#ax.set_zscale('log')
ax.tick_params(labelsize=12)
surf = ax.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount,
                       facecolors=colors, shade=False)
surf.set_facecolor((0,0,0,0))
plt.show()

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
