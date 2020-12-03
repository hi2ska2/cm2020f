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

nx = 61
ny = 15
N = nx*ny

#coefficients
mu = 0.15
Nc = 2.86e+25
kT = sc.k*300
ni = Nc*np.exp(-0.56*sc.e/kT)
Np = -1e+26
Nacc = 1e+26
q = sc.e 
e1 = 3.9*sc.epsilon_0
e2 = 11.7*sc.epsilon_0

interx1 = int((nx-1)/3.0)
interx2 = int(2*(nx-1)/3.0)
intery1 = int((ny-1)/7.0)
intery2 = int(6*(ny-1)/7.0)

VT = kT/q
Dn = VT*mu

Ndon = 1e26*np.ones(nx)
Ndon[interx1:interx2+1] = 0

dx = 30*1e-9/(nx-1)
deltax = 30*1e-9/(nx-1)
current = np.zeros(10)
voltage = np.zeros(10)

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

def init_phi(phi,n,GV,V_applied) :
	dphi_bc = np.zeros((nx*ny,nx*ny))
	#setting phi_bc and corresponding H
	H = np.zeros((nx*ny,nx*ny))
	phi_bc = np.zeros(nx*ny)
	#Boundary values on the gate line
	phi_bc[interx1:interx2+1] = 0.33374 + GV
	phi_bc[nx*(ny-1)+interx1:nx*(ny-1)+interx2+1] = 0.33374+GV

	np.fill_diagonal(H[0:nx*(ny-1)+interx2+1],1)
	np.fill_diagonal(H[0:nx*(ny-1)+interx1],0)
	np.fill_diagonal(H[0:interx2+1],1)
	np.fill_diagonal(H[0:interx1],0)

	#Boundary valeus on the silicon layer (left & right sides) except interlayer
	for i in range(intery2-intery1-1) :
		phi_bc[nx*int(intery1 + 1 + i)] = VT*np.arcsinh(Nacc/(ni*2.0))
		phi_bc[nx*int(intery1 + 1 + i + 1) -1] = V_applied + VT*np.arcsinh(Nacc/(ni*2.0))
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
		H[rd,rd] = -2*e1;H[rd,rd-nx] = 0.5*e1;H[rd,rd+nx]=0.5*e1; H[rd,rd-1]=e1;
		
	#interfaces
	lu = nx*intery1; ru = nx*intery1+nx-1;
	ld = nx*intery2; rd = nx*intery2+nx-1;

	H[lu,lu] = -(e1+e2); H[lu,lu-nx] = 0.5*e1; H[lu,lu+nx] = 0.5*e2; H[lu,lu+1] = (e1+e2)/2;
	H[ru,ru] = -(e1+e2); H[ru,ru-nx] = 0.5*e1; H[ru,ru+nx] = 0.5*e2; H[ru,ru-1] = (e1+e2)/2;
	H[ld,ld]= -(e1+e2);H[ld,ld-nx] = 0.5*e2;H[ld,ld+nx]=0.5*e1; H[ld,ld+1]=(e1+e2)/2;
	H[rd,rd]= -(e1+e2);H[rd,rd-nx] = 0.5*e2;H[rd,rd+nx]=0.5*e1; H[rd,rd-1]=(e1+e2)/2;

	#middle part
	for j in range(ny-2):
		for i in range(nx-2):
			centeru = nx*(j+1)+i+1
			if(j+1<intery1) : H[centeru,centeru] = -4*e1; H[centeru,centeru-1] = e1; H[centeru,centeru+1] = e1; H[centeru,centeru+nx]=e1;H[centeru,centeru-nx] = e1;
			elif(j+1==intery1) : H[centeru,centeru] = -2*(e1+e2); H[centeru,centeru-1] = (e1+e2)/2; H[centeru,centeru+1] = (e1+e2)/2; H[centeru,centeru+nx]=e2;H[centeru,centeru-nx] = e1;
			elif((j+1>intery1)&(j+1<intery2)) : 
				H[centeru,centeru] = -4*e2; H[centeru,centeru-1] = e2; H[centeru,centeru+1] = e2; H[centeru,centeru+nx]=e2; H[centeru,centeru-nx] = e2;
				if((i+1>=interx1)&(i+1<=interx2)):
					phi_bc[centeru] = dx*dx*q*n[centeru]
					dphi_bc[centeru,centeru] = 0
				else : 
					phi_bc[centeru] = dx*dx*q*(-Nacc+n[centeru])
					dphi_bc[centeru,centeru] = 0
			elif(j+1==intery2) : H[centeru,centeru] = -2*(e1+e2); H[centeru,centeru-1] = (e1+e2)/2; H[centeru,centeru+1] = (e1+e2)/2; H[centeru,centeru+nx]=e1;H[centeru,centeru-nx] = e2;
			else : H[centeru,centeru] = -4*e1; H[centeru,centeru-1] = e1; H[centeru,centeru+1] = e1; H[centeru,centeru+nx]=e1;H[centeru,centeru-nx] = e1;
	residue = np.matmul(H,phi) - phi_bc
	Rphi_phi = H -dphi_bc
	Rphi_n = np.zeros((N,N))
	#np.fill_diagonal(Rphi_n[0:nx*intery2],-dx*dx*q)
	#np.fill_diagonal(Rphi_n[0:nx*intery1],0)
	for i in range(nx*intery1,nx*intery2):
		x=i%nx
		y=i/nx
		if((x!=0)&(x!=nx-1)) :
			Rphi_n[i,i] = -dx*dx*q; 
	return residue, Rphi_phi, Rphi_n
'''
def density_cl(gate) :
	phi = np.full(nx*ny,0.33374)
	dphi_bc = np.zeros((nx*ny,nx*ny))

	#setting phi_bc and corresponding H
	H = np.zeros((nx*ny,nx*ny))
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
					dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*(q/kT)

				else : 
					phi_bc[centeru] = dx*dx*q*(Nacc+2*ni*np.sinh(q*phi[centeru]/kT))
					dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*(q/kT)
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

	#centerphi = np.zeros(39)
	#while abs(tmp-min_phi)>1e-5 :
	for i in range(20):
		print("itr : ",itr+1,"abs(diff) : ",tmp-min_phi)
		min_phi = np.min(phi)
		residue = np.matmul(H,phi) - phi_bc
		Jaco = H -dphi_bc
		phi = phi - slin.solve(Jaco,residue)
		tmp = np.min(phi)
		#centerphi[itr] = abs(tmp-min_phi)
		centerphi = abs(tmp-min_phi)
		for j in range(ny-2):
			for i in range(nx-2):
				centeru = nx*(j+1)+i+1
				if((j+1>intery1)&(j+1<intery2)) : 
					if((i+1>=interx1)&(i+1<=interx2)):
						phi_bc[centeru] = dx*dx*q*2*ni*np.sinh(q*phi[centeru]/kT)
						dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*(q/kT)
					else : 
						phi_bc[centeru] = dx*dx*q*(Nacc+2*ni*np.sinh(q*phi[centeru]/kT))
						dphi_bc[centeru,centeru] = dx*dx*q*2*ni*np.cosh(q*phi[centeru]/kT)*(q/kT)
		itr = itr+1			
	return phi
'''

def init_e(phi,n):
	residue = copy.deepcopy(n)
	Rn_n=np.zeros((N,N))
	for i in range(N):
		Rn_n[i,i] = 1
	Rn_phi=np.zeros((N,N))
	#silicon layer
	for i in range(nx*(intery1+1),nx*intery2):
		x = i%nx; y=int(i/nx);
		if x==0:
			residue[i] = n[i] - Nacc
		elif (x==nx-1) :
			residue[i] = n[i] - Nacc
		else:
			residue[i] = n[i+1]*Bn((phi[i+1]-phi[i])/VT) - n[i]*Bn((phi[i]-phi[i+1])/VT) - n[i]*Bn((phi[i]-phi[i-1])/VT) + n[i-1]*Bn((phi[i-1]-phi[i])/VT)
			Rn_n[i,i] = -Bn((phi[i]-phi[i+1])/VT) - (Bn((phi[i]-phi[i-1])/VT))
			Rn_n[i,i+1] = Bn((phi[i+1]-phi[i])/VT)
			Rn_n[i,i-1] = Bn((phi[i-1]-phi[i])/VT)	
			
			Rn_phi[i,i] = -n[i+1]*dBn((phi[i+1]-phi[i])/VT)/VT-n[i]*dBn((phi[i]-phi[i+1])/VT)/VT - (n[i]*dBn((phi[i]-phi[i-1])/VT)/VT+n[i-1]*dBn((phi[i-1]-phi[i])/VT)/VT)
			Rn_phi[i,i+1] = n[i+1]*dBn((phi[i+1]-phi[i])/VT)/VT+n[i]*dBn((phi[i]-phi[i+1])/VT)/VT 
			Rn_phi[i,i-1] = n[i]*dBn((phi[i]-phi[i-1])/VT)/VT+n[i-1]*dBn((phi[i-1]-phi[i])/VT)/VT
	return residue, Rn_n, Rn_phi

def merge_jaco(Rphi_phi,Rphi_n,Rn_n,Rn_phi) :
	Jaco = np.zeros((2*N,2*N))
	Jaco[::2,::2] = Rphi_phi
	Jaco[::2,1::2] = Rphi_n
	Jaco[1::2,::2] = Rn_phi
	Jaco[1::2,1::2] = Rn_n
	return Jaco

def scaling_update(Jaco, residue):
	Cvector = np.zeros(2*N)
	Cvector[0::2] = VT 
	Cvector[1::2] = np.max(np.abs(Ndon))
	Cmatrix = np.diag(Cvector)
	Jaco_scaled = np.matmul(Jaco,Cmatrix)
	Rvector = 1./np.sum(np.abs(Jaco_scaled),axis=1)
	Rmatrix = np.diag(Rvector)
	Jaco_scaled = np.matmul(Rmatrix,Jaco_scaled)
	res_scaled = np.matmul(Rmatrix,residue)
	#plt.matshow(Jaco_scaled[::2,::2])
	#plt.show()
	update_scaled = slin.solve(Jaco_scaled,-res_scaled)
	update = np.matmul(Cmatrix,update_scaled)
	return update

def Jn(i,j,n,phi):
	return n[i]*Bn((phi[i]-phi[j])/VT)-n[j]*Bn((phi[j]-phi[i])/VT)

def current(n, phi):
	curry=0
	for i in range(intery1,intery2-2):
		curry = curry+np.abs(q*Dn*Jn((i+1+1)*nx-1,(i+1+1)*nx-2 ,n,phi)/dx)
	return curry

def DD(phi,n, gate) :
	phi = np.zeros(N)
	GV=0
	for gate in range(1):
		I = np.zeros(12)
		n = np.ones(N)
		phi = np.zeros(N)
		GV = 0.1*gate
		GV = 0.5
		V_applied = 0
		for applied in range(1):
			V_applied = 0.1*applied
			V_applied = 0.5
			print("VG,Vapp = (",GV,V_applied,")")
			for i in range(5):
				print("Jaco step : ",i)
				Jaco = np.zeros((2*N,2*N))
				residue_e = np.zeros(N)
				residue_phi = np.zeros(N)
				residue = np.zeros(2*N)
				residue_phi, Rphi_phi, Rphi_n = init_phi(phi,n,GV,V_applied)
				residue_e, Rn_n, Rn_phi = init_e(phi,n)
				Jaco = merge_jaco(Rphi_phi,Rphi_n,Rn_n,Rn_phi)
				residue[::2] = residue_phi; residue[1::2] = residue_e;
				delta = scaling_update(Jaco,residue)
				n = n + delta[1::2]
				phi = phi + delta[::2]
			I[applied] = current(n,phi)
		V = np.linspace(0,1.1,num=12)
		plt.plot(V,I,c=plt.cm.viridis(gate/11.0),linestyle='dashed',label='$V_{gate}=$'+str(round(GV,3))+'V',marker='o',lw=1,ms=10)
		plt.xlabel('$V_D$',fontsize=18)
		plt.ylabel('I',fontsize=18)
		plt.legend(fontsize=18)
		plt.tick_params(labelsize=18)
	return phi, n
		
gate = 0
n = np.ones(N)
phi = np.zeros(N)
phi, n = DD(phi,n,gate)
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

Z = np.reshape(n,(-1,nx))
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
ax.set_zlabel(r"$n(x,y)$",fontsize=12,rotation=90)
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
