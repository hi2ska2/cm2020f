import matplotlib.pyplot as plt
import numpy as np
import copy
import scipy.constants as sc
import scipy.linalg as slin
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import scipy.integrate as inte
import fractions
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.tri as mtri
import matplotlib.ticker as mticker

x, y = np.loadtxt("MERGED.vertex",usecols=(0,1),unpack=True)
fst,scd,thd = np.loadtxt("MERGED.element",usecols=(0,1,2),unpack=True)
NofTri = np.size(fst)
N = np.size(x)
idx = np.zeros((NofTri,3))
idx[:,0] = fst-1; idx[:,1] = scd-1; idx[:,2] = thd-1;
idx.astype(int)

f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.4e' % x))
fmt = mticker.FuncFormatter(g)


def oisim(i1,i2,i3):
	A1 = np.zeros((2,2))
	A1[0,0] = x[i2] - x[i1];  A1[0,1] = y[i2] - y[i1]; 
	A1[1,0] = x[i3] - x[i2];  A1[1,1] = y[i3] - y[i2]; 
	B1 = np.zeros(2)
	B1[0] = x[i2]**2-x[i1]**2+y[i2]**2-y[i1]**2
	B1[1] = x[i3]**2-x[i2]**2+y[i3]**2-y[i2]**2
	B1 = B1/2.0
	oisim_x,oisim_y = slin.solve(A1,B1)
	return oisim_x, oisim_y

def makeLA(idx,x,y) :
	L = np.zeros((N,N))
	A = np.zeros((N,N))
	#the length of edges are calculated here
	for i in range(NofTri):
		i1 = int(idx[i,0]); i2 = int(idx[i,1]); i3 = int(idx[i,2]);
		o_x, o_y = oisim(i1,i2,i3)
		#plt.scatter(o_x,o_y,c='b')
		R = np.sqrt((o_x-x[i1])**2+(o_y-y[i1])**2)
		L[i1,i2] = np.sqrt((x[i1]-x[i2])**2+(y[i1]-y[i2])**2)
		L[i2,i3] = np.sqrt((x[i2]-x[i3])**2+(y[i2]-y[i3])**2)
		L[i1,i3] = np.sqrt((x[i1]-x[i3])**2+(y[i1]-y[i3])**2)
		L[i2,i1] = L[i1,i2]; L[i3,i2] = L[i2,i3]; L[i3,i1] = L[i1,i3];
		'''
		if(i<10) : 
			plt.text((x[i1]+x[i2])/2,(y[i1]+y[i2])/2,'L={}'.format(fmt(L[i1,i2])),fontsize=18)
			plt.text((x[i2]+x[i3])/2,(y[i2]+y[i3])/2,'L={}'.format(fmt(L[i2,i3])),fontsize=18)
			plt.text((x[i3]+x[i1])/2,(y[i3]+y[i1])/2,'L={}'.format(fmt(L[i1,i3])),fontsize=18)
		'''
		if R**2 - (L[i1,i2]/2.0)**2 > 0 :
			A[i1,i2] = np.sqrt(R**2 - (L[i1,i2]/2.0)**2)
		if R**2 - (L[i2,i3]/2.0)**2 > 0 :
			A[i2,i3] = np.sqrt(R**2 - (L[i2,i3]/2.0)**2)
		if R**2 - (L[i1,i3]/2.0)**2 > 0 :
			A[i1,i3] = np.sqrt(R**2 - (L[i1,i3]/2.0)**2)
		A[i2,i1] = A[i1,i2]; A[i3,i2] = A[i2,i3]; A[i3,i1] = A[i1,i3];
	return L,A

def makeH (L,A) : 
	C = np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			if L[i,j] >1e-11 : C[i,j] = A[i,j]/L[i,j]
	H = np.zeros((N,N))
	for i in range(N):
		for j in range(N):
			if C[i,j]!=0 :
				H[i,i] = H[i,i] -C[i,j]
				H[i,j] = C[i,j]
	return H

triang = mtri.Triangulation(x,y,triangles = idx)
#triang = mtri.Triangulation(x,y)
#plt.show()

'''
for i in range(N):
	if(abs(y[i]) < 1e-9) : 
		print('x = ',x[i],' y = ',y[i],' idx = ',i)

for i in range(N):
	if(abs(y[i]) == 1.85e-7) : 
		print('x = ',x[i],' y = ',y[i],' idx = ',i)
'''

for i in range(N):
	if (x[i]<3.8e-8)&(x[i] > 3.6e-8) : 
		print('x = ',x[i],' y = ',y[i],' idx = ',i)


#L is the length bet edges, A is the Boronoi Area bet edges, C = A/L, H is the coefficients to solve Laplace eqn.
L ,A = makeLA(idx,x,y)

#plt.plot(V,I,c=plt.cm.viridis(gate/11.0),linestyle='dashed',label='$V_{gate}=$'+str(round(GV,3))+'V',marker='o',lw=1,ms=10)
#plt.legend(fontsize=18)

H = makeH(L,A)

print("We have ",N," points. You have to set Phi to initiate eqns.")
index_0 = int(input("Phi = 0, Enter index : "))
index_1 = int(input("Phi = 1, Enter index : "))

H[index_0,:] = 0; H[index_0,index_0] = 1;
H[index_1,:] = 0; H[index_1,index_1] = 1;

#for i in range(N):
#	if H[1001,i]!=0 : print(H[1001,i])

plt.triplot(triang,marker='o',ms=3,c='k',mfc='r',zorder=1)
t = plt.annotate(r'$\phi=0$',xy=(x[index_0],y[index_0]),fontsize=15,xytext=(x[index_0]-1e-7,y[index_0]),arrowprops=dict(facecolor='b',shrink=0.05),zorder=2)
t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='white'))

plt.scatter(x[index_0],y[index_0],c='b',s=80,zorder=3)
plt.annotate(r'$\phi=1$',xy=(x[index_1],y[index_1]),fontsize=15,xytext=(x[index_1]+1e-7,y[index_1]),arrowprops=dict(facecolor='r',shrink=0.05),zorder=2)
plt.scatter(x[index_1],y[index_1],c='r',s=80,zorder=3)
plt.tick_params(labelsize=18)
plt.xlabel('$x$',fontsize=18)
plt.ylabel('$y$',fontsize=18)
plt.show()
boundary = np.zeros(N)
boundary[index_0] = 0
boundary[index_1] = 1

phi = slin.solve(H,boundary)

Z = np.reshape(phi,(-1,N))
X,Y = np.meshgrid(x,y)

#Z = np.log(Z)
norm = plt.Normalize(Z.min(), Z.max())
colors = cm.viridis(norm(Z))
rcount, ccount, _ = colors.shape

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.zaxis.set_rotate_label(False)
ax.text(x[index_0],y[index_0],0,'(x,y,z) =({},{},0)'.format(fmt(x[index_0]),fmt(y[index_0]),fontsize=18))
ax.text(x[index_1],y[index_1],1,'(x,y,z) =({},{},1)'.format(fmt(x[index_1]),fmt(y[index_1]),fontsize=18))
ax.set_xlabel("x",fontsize=12)
ax.set_ylabel("y",fontsize=12)
ax.set_zlabel(r"$\phi(x,y)$",fontsize=12,rotation=90)
#ax.set_zscale('symlog')
ax.tick_params(labelsize=12)
surf = ax.plot_trisurf(triang,phi,cmap=plt.cm.jet)
#surf = ax.plot_trisurf(x,y,phi,cmap=plt.cm.jet)
#surf = ax.plot_surface(X, Y, Z, rcount=rcount, ccount=ccount,
#                       facecolors=colors, shade=False)
surf.set_facecolor((0,0,0,0))
plt.show()


