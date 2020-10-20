import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions

#np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})

Nc = 2.86*1e+25
#kT = 25.85*1e-3
kT = sc.k*300
ni = Nc*np.exp(-0.56*sc.e/kT)
N = 1000
Nacc = 1e+24
q = sc.e 
e1 = 3.9*sc.epsilon_0
e2 = 11.7*sc.epsilon_0
deltax = ((6.6*1e-9)/(N-1))
Lx = 1e-7
Ly = 1e-7
mzz = 0.19*sc.m_e

inter2 = int((N-1)*(5.8/6.6))
inter1 = int((N-1)*(0.8/6.6))
thermal = kT/q

itr = 20 
centerphi = np.zeros(itr)
def density_phi_cl(gate):
    p_density=np.zeros(N)
    e_density=np.zeros(N)
    density=np.zeros(N)
    phi = np.full(N,0.33374)
    #for i in range(itr) :
    GV = (gate)/10.0
    centerphi = 1000
    diff = 1000
    while diff>1e-5 :
        centerphi = phi[int(N/2)]
        e_density = ni*np.exp(q*phi/kT)
        e_density[:inter1] = 0
        e_density[inter2+1:] = 0

        #p_density = ni*np.exp(-q*phi/kT)
        #p_density[:inter1] = 0
        #p_density[inter2+1:] = 0
 
        density = 2*ni*np.sinh(q*phi/kT)
        density[:inter1] = 0
        density[inter2+1:] = 0

        H = np.zeros((N,N))
        np.fill_diagonal(H,-2*e1)
        np.fill_diagonal(H[1:],e1)
        np.fill_diagonal(H[:,1:],e1)

        np.fill_diagonal(H[0:inter2+1],-2*e2)
        np.fill_diagonal(H[1:inter2+1],e2)
        np.fill_diagonal(H[:,1:inter2+1],e2)
        H[inter2][inter2] = (-e2-e1)

        np.fill_diagonal(H[0:inter1+1],-2*e1)
        np.fill_diagonal(H[1:inter1+1],e1)
        np.fill_diagonal(H[:,1:inter1+1],e1)
        H[inter1][inter1] = (-e2-e1)

        H[0][0] = 1
        H[0][1] = 0
        H[N-1][N-1] = 1
        H[N-1][N-2] = 0

        B = np.zeros(N)
        B[inter1] = q*(Nacc+density[inter1])/2.0
        B[inter1+1:inter2] = q*(Nacc+density[inter1+1:inter2])
        B[inter2] = q*(Nacc+density[inter2])/2.0
        B = B*deltax*deltax
        B[0] = 0.33374 + GV
        B[N-1] = 0.33374 + GV

        res = np.zeros(N)
        res = np.matmul(H,phi) - B

        der_B = np.zeros((N,N))
        np.fill_diagonal(der_B[inter1:inter2,inter1:inter2],2*q*(q/kT)*ni*np.cosh(q*phi[inter1:inter2]/kT))
        der_B[inter1,inter1] = (q/kT)*q*ni*np.cosh(q*phi[inter1]/kT)
        der_B[inter2,inter2] = (q/kT)*q*ni*np.cosh(q*phi[inter2]/kT)

        der_B = der_B*deltax*deltax

        Jaco = np.zeros((N,N))
        Jaco = H - der_B
        
        phi=phi + slin.solve(Jaco,-res)

        diff = abs(phi[int(N/2)]-centerphi)
    x = np.linspace(0,6.6*1e-9,num=N,endpoint=True)
    #plt.plot(x, density,label=r'$V_{gate}=$'+str(GV)+'V',c=plt.cm.cool(gate/11),marker='o',lw=0,ms=10)
    return density,phi

def density_cooking(phi) :
	bandgap_si = 1.11
	bandgap_o = 8.9

	H = np.zeros((N-2,N-2))
	np.fill_diagonal(H,sc.hbar*sc.hbar/(mzz*deltax*deltax)-q*phi+q*bandgap_o/2.0)
	np.fill_diagonal(H[:inter2-1],sc.hbar*sc.hbar/(mzz*deltax*deltax)-q*phi+q*bandgap_si/2.0)
	np.fill_diagonal(H[:inter1-1],sc.hbar*sc.hbar/(mzz*deltax*deltax)-q*phi+q*bandgap_o/2.0)
	np.fill_diagonal(H[1:],-sc.hbar*sc.hbar/(2*mzz*deltax*deltax))
	np.fill_diagonal(H[:,1:],-sc.hbar*sc.hbar/(2*mzz*deltax*deltax))

	E,psi = slin.eigh(H)
	psi = psi/np.sqrt(deltax)
	psi = np.concatenate((np.zeros((1,N-2)),psi),axis=0)
	psi = np.concatenate((psi,np.zeros((1,N-2))),axis=0)

	#intagrating ground state -> should be 1.
	z = np.linspace(0,6.6*1e-9,N,endpoint=True)
	print(inte.simps(psi[:,0]**2,z,deltax))

	density = 2*np.sum((psi**2)/((1+np.exp(E/kT))*Lx*Ly),axis=1)
	return density

def density_phi_qm(gate,density,phi):
    GV = (gate)/10.0
    diff = 1000
    centerphi = 1000
    while diff>1e-15:
        centerphi = phi[int(N/2)]
        H = np.zeros((N,N))
        np.fill_diagonal(H,-2*e1)
        np.fill_diagonal(H[1:],e1)
        np.fill_diagonal(H[:,1:],e1)

        np.fill_diagonal(H[0:inter2+1],-2*e2)
        np.fill_diagonal(H[1:inter2+1],e2)
        np.fill_diagonal(H[:,1:inter2+1],e2)
        H[inter2][inter2] = (-e2-e1)

        np.fill_diagonal(H[0:inter1+1],-2*e1)
        np.fill_diagonal(H[1:inter1+1],e1)
        np.fill_diagonal(H[:,1:inter1+1],e1)
        H[inter1][inter1] = (-e2-e1)

        H[0][0] = 1
        H[0][1] = 0
        H[N-1][N-1] = 1
        H[N-1][N-2] = 0

        B = np.zeros(N)
        B[inter1] = q*(Nacc+density[inter1])/2.0
        B[inter1+1:inter2] = q*(Nacc+density[inter1+1:inter2])
        B[inter2] = q*(Nacc+density[inter2])/2.0
        B = B*deltax*deltax
        B[0] = 0.33374 + GV
        B[N-1] = 0.33374 + GV
        phi = slin.solve(H,B)
        diff = abs(phi[int(N/2)]-centerphi)
        print(diff)
        density = density_cooking(phi)
    x = np.linspace(0,6.6*1e-9,num=N,endpoint=True)
    plt.plot(x, density,label=r'$V_{gate}=$'+str(GV)+'V',c=plt.cm.plasma(gate/11),marker='o',lw=0,ms=10)
    return density,phi

for gate in range(11) :
    density, phi = density_phi_cl(gate)
    q_density = density_cooking(phi)
    q_density, phi = density_phi_qm(gate,q_density,phi)

    plt.xlabel("z",fontsize = 18)
    plt.ylabel(r"$n(z)$",fontsize = 18)
    plt.legend(fontsize=18)
    #plt.yscale('symlog')
    plt.tick_params(labelsize=18)
    plt.tight_layout()
    

plt.show()

