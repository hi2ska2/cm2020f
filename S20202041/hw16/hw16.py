import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions

kT = sc.k*300
ni = 1.075*1e+16
N = 601
q = sc.e 
e2 = 11.7
deltax = ((60*1e-9)/(N-1))
mu = 1.5*1e-1

inter1 = int((N-1)/6.0)
inter2 = int(5*(N-1)/6.0)
thermal = kT/q
VT = kT/q

Ndon = 5e23*np.ones(N)
Ndon[inter1:inter2+1] = 2e21
current = np.zeros(10)
voltage = np.zeros(10)


def density_phi_cl(gate):
    e_density=np.zeros(N)
    e_density[0] = Ndon[0]; e_density[N-1] = Ndon[N-1]
    density=np.zeros(N)
    phi = kT*np.log(Ndon/ni)/q
    #for i in range(itr) :
    GV = (gate)/10.0
    centerphi = 1000
    diff = 1000
    while diff>1e-15 :
        centerphi = phi[int(N/2)]
        e_density = ni*np.exp(q*phi/kT)
        e_density[0] = Ndon[0]; e_density[N-1] = Ndon[N-1]

        density = 2*ni*np.sinh(q*phi/kT)

        H = np.zeros((N,N))
        np.fill_diagonal(H,-2*e2*sc.epsilon_0)
        np.fill_diagonal(H[1:],e2*sc.epsilon_0)
        np.fill_diagonal(H[:,1:],e2*sc.epsilon_0)

        H[0][0] = 1
        H[0][1] = 0
        H[N-1][N-1] = 1
        H[N-1][N-2] = 0

        B = np.zeros(N)
        #B[:] = q*(-Ndon+density)/2.0
        B[:] = q*(-Ndon+density)
        B[inter1:inter2+1] = q*(-Ndon[inter1:inter2+1]+density[inter1:inter2+1])
        B = B*deltax*deltax
        B[0] = phi[0] ; B[N-1] = phi[N-1];

        res = np.zeros(N)
        res = np.matmul(H,phi) - B

        der_B = np.zeros((N,N))
        np.fill_diagonal(der_B[:,:],2*deltax*deltax*q*(q/kT)*ni*np.cosh(q*phi[:]/kT))
        #der_B[inter,inter] = (q/kT)*q*ni*np.cosh(q*phi[inter]/kT)
        #np.fill_diagonal(der_B[:inter,:inter],2*q*(q/kT)*ni*np.cosh(q*phi[:inter]/kT))

        Jaco = np.zeros((N,N))
        Jaco = H - der_B
        
        phi=phi + slin.solve(Jaco,-res)

        diff = abs(phi[int(N/2)]-centerphi)
    return phi, e_density

def DD(phi,e_density):
    for bias in range(10):
	    V_app = 0.05*bias
	    centerphi = 1000
	    diff = 1000
	    while(diff>1e-10):
		    Jaco = np.zeros((2*N,2*N))
		    residue_e = np.zeros(N)
		    residue_phi=np.zeros(N)
		    residue = np.zeros(2*N)

		    der_phi = np.zeros(N)
		    der_e = np.zeros(N); avg_e = np.zeros(N)
		 

		    centerphi = phi[int(N/2)]
		    der_phi[:N-1] = (phi[1:]-phi[:N-1])/deltax
		    der_e[:N-1] = (e_density[1:]-e_density[:N-1])/deltax
		    avg_e[:N-1] = (e_density[1:]+e_density[:N-1])/2.0

		    Jn = avg_e*der_phi-VT*der_e
		    dJn = der_phi/2.0+VT/deltax
		    dJn_ = der_phi/2.0-VT/deltax

		    residue_phi[1:N-1] = e2*phi[2:N] - 2*e2*phi[1:N-1]+e2*phi[0:N-2]+deltax*deltax*q*(Ndon[1:N-1]-e_density[1:N-1])/sc.epsilon_0
		    residue_e[1:] = Jn[1:] - Jn[:N-1]

		    residue[::2] = residue_phi
		    residue[1::2] = residue_e
		    #residue[::2] = e2*phi[2:N] - 2*e2*phi[1:N-1]+e2*phi[0:N-2]+deltax*deltax*q*(Ndon[1:N-1]-ni)
		    #residue[1::2] = Jn[1:] - Jn[:N-1]

		    #allocatong derivative of R_phi by phi
		    np.fill_diagonal(Jaco[::2,::2],-2*e2)
		    np.fill_diagonal(Jaco[2::2,::2],e2)
		    np.fill_diagonal(Jaco[::2,2::2],e2)
		    Jaco[2+2*(N-2),0+2*(N-2)] = 0

		    #allocatong derivative of R_phi by n
		    np.fill_diagonal(Jaco[2::2,3::2],-deltax*deltax*q/sc.epsilon_0)
		    Jaco[0+2*(N-1),1+2*(N-1)] = 0; 	    
		    #allocatong derivative of R_n by n
		    np.fill_diagonal(Jaco[3::2,3::2],dJn[1:]-dJn_[:N-1])
		    np.fill_diagonal(Jaco[1::2,3::2],dJn_[0:])
		    np.fill_diagonal(Jaco[3::2,1::2],-dJn[:N-1])
		    #allocatong derivative of R_n by phi
		    np.fill_diagonal(Jaco[3::2,4::2],avg_e[1:]/deltax)
		    np.fill_diagonal(Jaco[3::2,2::2],-avg_e[1:]/deltax-avg_e[:N-1]/deltax)
		    np.fill_diagonal(Jaco[3::2,::2],avg_e[:N-1]/deltax)
		    Jaco[1+2*(N-1),0+2*(N-1)] = 0; Jaco[3+2*(N-2),0+2*(N-2)] = 0
		    #boundary condition by phi   
		    residue[0] = phi[0]-VT*np.log(Ndon[0]/ni)
		    residue[2*N-2] =phi[N-1]-VT*np.log(Ndon[N-1]/ni)-V_app
		    #Jaco[0,:] = 0; 
		    Jaco[0,0] = 1
		    #Jaco[2*N-2,:] = 0; 
		    Jaco[2*N-2,2*N-2] = 1  
		    #boundary condition by n
		    residue[1] = e_density[0]-Ndon[0]
		    residue[2*N-1] = e_density[N-1]-Ndon[N-1]
		    Jaco[1,:]=0; Jaco[1,1] = 1 ;
		    Jaco[2*N-1,:] = 0 ; Jaco[2*N-1,2*N-1] = 1

		    delta = slin.solve(Jaco,-residue)
		    e_density = e_density + delta[1::2]
		    phi = phi + delta[::2]
		    diff = abs(phi[int(N/2)]-centerphi)
		    print(diff)
	    plt.plot(x, e_density,label='$V_{app}=$'+str(round(V_app,3))+'V',c=plt.cm.plasma(bias/10),marker='o',lw=0,ms=10)
	    current[bias] = Jn[N-2]*q*mu/1e-4
	    voltage[bias] = V_app
    return phi,e_density

gate = 0
x = np.linspace(0,60*1e-9,num=N,endpoint=True)
phi, e_density = density_phi_cl(gate)

old_e = e_density; old_phi = phi

phi, e_density = DD(phi,e_density)
#plt.plot(x, e_density,label=r'$\phi_{DD}$',c='sienna',marker='s',lw=0,ms=10)
#plt.plot(x, e_density,label=r'$n_{DD}$',c='k',marker='s',lw=0,ms=10,mfc='none')
plt.xlabel("z",fontsize = 18)
#plt.ylabel(r"$\phi(z)$",fontsize = 18)
plt.ylabel(r"$n(z)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('symlog')
plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()

plt.plot(voltage, current,label=r'$n$',c='r',marker='o',lw=0,ms=10,mfc='none')
plt.xlabel(r"$V_{applied}$",fontsize = 18)
plt.ylabel(r"$I$",fontsize = 18)
#plt.legend(fontsize=18)
#plt.yscale('symlog')
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
'''
plt.plot(x, old_phi,label=r'$\phi$',c='r',marker='o',lw=0,ms=12,mfc='none')
plt.plot(x, phi,label=r'$\phi_{DD}$',c='k',marker='s',lw=0,ms=10,mfc='none')
plt.xlabel("z",fontsize = 18)
plt.ylabel(r"$\phi(z)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('symlog')
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
'''
