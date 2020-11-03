import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions

Nc = 2.86*1e+25
kT = sc.k*300
ni = 1.075*1e+16
N = 601
#N = 51
q = sc.e 
e2 = 11.7#*sc.epsilon_0
deltax = ((600*1e-9)/(N-1))

inter1 = int((N-1)/6.0)
inter2 = int(5*(N-1)/6.0)
thermal = kT/q
Ndon = 5e23*np.ones(N)
Ndon[inter1:inter2+1] = 2e21


def density_phi_cl(gate):
    e_density=np.zeros(N)
    e_density[0] = Ndon[0]; e_density[N-1] = Ndon[N-1]
    density=np.zeros(N)
    phi = kT*np.log(Ndon/ni)/q
    #for i in range(itr) :
    GV = (gate)/10.0
    centerphi = 1000
    diff = 1000
    while diff>1e-8 :
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
        B[:inter1] = q*(-Ndon[:inter1]+density[:inter1])/2.0
        B[inter1:inter2+1] = q*(-Ndon[inter1:inter2+1]+density[inter1:inter2+1])/2.0
        B[inter2+1:] = q*(-Ndon[inter2+1:]+density[inter2+1:])/2.0
        B = B*deltax*deltax
        B[0] = phi[0] ; B[N-1] = phi[N-1];

        res = np.zeros(N)
        res = np.matmul(H,phi) - B

        der_B = np.zeros((N,N))
        np.fill_diagonal(der_B[:,:],2*q*(q/kT)*ni*np.cosh(q*phi[:]/kT))
        #der_B[inter,inter] = (q/kT)*q*ni*np.cosh(q*phi[inter]/kT)
        #np.fill_diagonal(der_B[:inter,:inter],2*q*(q/kT)*ni*np.cosh(q*phi[:inter]/kT))

        der_B = der_B*deltax*deltax

        Jaco = np.zeros((N,N))
        Jaco = H - der_B
        
        phi=phi + slin.solve(Jaco,-res)

        diff = abs(phi[int(N/2)]-centerphi)
    return phi, e_density

def DD(phi,e_density):
    VT = kT/q
    Jaco = np.zeros((2*N,2*N))
    Jaco_phi_phi = np.zeros((N,N))
    Jaco_phi_n = np.zeros((N,N))
    Jaco_n_phi = np.zeros((N,N))
    Jaco_n_n = np.zeros((N,N))

    residue_e = np.zeros(N)
    residue_phi=np.zeros(N)
    residue = np.zeros(2*N)

    der_phi = np.zeros(N)
    der_e = np.zeros(N); avg_e = np.zeros(N)
    centerphi = 1000
    diff = 1000
    #while(diff>1e-8):
    for i in range(10):
	    centerphi = phi[int(N/2)]
	    der_phi[:N-1] = (phi[1:]-phi[:N-1])/deltax
	    der_e[:N-1] = (e_density[1:]-e_density[:N-1])/deltax
	    avg_e[:N-1] = (e_density[1:]+e_density[:N-1])/2.0

	    Jn = avg_e*der_phi-VT*der_e
	    dJn = der_phi/2.0+VT/deltax
	    dJn_ = der_phi/2.0-VT/deltax

	    phi0=kT*np.log(5e23/ni)/q
	    H = np.zeros((N,N))

	    np.fill_diagonal(H,-2.*e2)
	    np.fill_diagonal(H[1:],e2)
	    np.fill_diagonal(H[:,1:],e2)

	    Rphi = np.matmul(H,phi)+q/sc.epsilon_0*(Ndon-e_density)*deltax**2
	    Rphi[0] = phi[0]-phi0; Rphi[N-1] = phi[N-1]-phi0
	    
	    dphi = np.zeros(N); dphi[:N-1] = (phi[1:]-phi[:N-1])/deltax
	    dn = np.zeros(N); dn[:N-1] = (e_density[1:]-e_density[:N-1])/deltax
	    nav = np.zeros(N); nav[:N-1] = (e_density[1:]+e_density[:N-1])/2.

	    Jn = nav*dphi-VT*dn
	    Rn = np.zeros(N); Rn[1:] = Jn[1:]-Jn[:N-1]
	    Rn[0] = e_density[0]-5e23; Rn[N-1] = e_density[N-1]-5e23

	    R = np.zeros(2*N)
	    R[0::2] = Rphi
	    R[1::2] = Rn


	    residue_phi[1:N-1] = e2*phi[2:N] - 2*e2*phi[1:N-1]+e2*phi[0:N-2]+deltax*deltax*q*(Ndon[1:N-1]-e_density[1:N-1])/sc.epsilon_0
	    residue_e[1:] = Jn[1:] - Jn[:N-1]

	    residue[::2] = residue_phi
	    residue[1::2] = residue_e
	    #residue[::2] = e2*phi[2:N] - 2*e2*phi[1:N-1]+e2*phi[0:N-2]+deltax*deltax*q*(Ndon[1:N-1]-ni)
	    #residue[1::2] = Jn[1:] - Jn[:N-1]

	    #allocatong derivative of R_phi by phi
	    np.fill_diagonal(Jaco_phi_phi,-2*e2)
	    np.fill_diagonal(Jaco_phi_phi[:,1:],e2)
	    np.fill_diagonal(Jaco_phi_phi[1:,:],e2)
	    
	    Jaco_phi_phi[N-1,N-1] = 1 ; Jaco_phi_phi[N-1,N-2] = 0;
	    Jaco_phi_phi[0,1] = 0 ; Jaco_phi_phi[0,0] = 1;
	    #allocatong derivative of R_phi by n
	    np.fill_diagonal(Jaco_phi_n,-deltax*deltax*q/sc.epsilon_0)
	    Jaco_phi_n[N-1,N-1] = 0 ; Jaco_phi_n[N-1,N-2] = 0;
	    Jaco_phi_n[0,1] = 0 ; Jaco_phi_n[0,0] = 0;
	    #allocatong derivative of R_n by n
	    np.fill_diagonal(Jaco_n_n[1:,1:],dJn[1:]-dJn_[:N-1])
	    np.fill_diagonal(Jaco_n_n[:,1:],dJn_)
	    np.fill_diagonal(Jaco_n_n[1:,:],-dJn)
	    Jaco_n_n[N-1,N-1] = 1 ; Jaco_n_n[N-1,N-2] = 0;
	    Jaco_n_n[0,1] = 0 ; Jaco_n_n[0,0] = 1;

	    #allocatong derivative of R_n by phi
	    np.fill_diagonal(Jaco_n_phi[1:],avg_e/deltax)
	    np.fill_diagonal(Jaco_n_phi[1:,1:],-avg_e[1:]/deltax-avg_e[:N-1]/deltax)
	    np.fill_diagonal(Jaco_n_phi[:,1:],avg_e/deltax)
	    Jaco_n_phi[N-1,N-1] = 0 ; Jaco_n_phi[N-1,N-2] = 0;
	    Jaco_n_phi[0,1] = 0 ; Jaco_n_phi[0,0] = 0;
	    #boundary condition by phi   
	    residue[0] = phi[0]-kT*np.log(Ndon[0]/ni)/q
	    residue[2*N-2] =phi[N-1]-kT*np.log(Ndon[N-1]/ni)/q
	    #boundary condition by n
	    residue[1] = e_density[0]-Ndon[0]
	    residue[2*N-1] = e_density[N-1]-Ndon[N-1]

	    Jaco[0::2,0::2] = Jaco_phi_phi
	    Jaco[0::2,1::2] = Jaco_phi_n
	    Jaco[1::2,0::2] = Jaco_phi_n
	    Jaco[1::2,1::2] = Jaco_n_n
 
	    delta = slin.solve(Jaco,-residue)
	    e_density = e_density + delta[1::2]
	    phi = phi + delta[::2]
	    diff = abs(phi[int(N/2)]-centerphi)
	    print(diff)


    '''
    while(diff>1e-8):
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
	    residue[0] = phi[0]-kT*np.log(Ndon[0]/ni)/q
	    residue[2*N-2] =phi[N-1]-kT*np.log(Ndon[N-1]/ni)/q
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
    '''
    print(Jaco[3,0:5])
    #plt.imshow(Jaco)
    #plt.show()
    return phi,e_density

gate = 0
x = np.linspace(0,600*1e-9,num=N,endpoint=True)
phi, e_density = density_phi_cl(gate)
plt.plot(x, e_density,label=r'$\phi$',c='r',marker='o',lw=0,ms=10)

old_e = e_density; old_phi = phi

phi, e_density = DD(phi,e_density)
#plt.plot(x, e_density,label=r'$\phi_{DD}$',c='sienna',marker='s',lw=0,ms=10)
plt.plot(x, e_density,label=r'$n_{density}$',c='b',marker='o',lw=0,ms=10)
plt.xlabel("z",fontsize = 18)
plt.ylabel(r"$\phi(z)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('symlog')
plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()




