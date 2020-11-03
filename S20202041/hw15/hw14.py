import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions

Nc = 2.86*1e+25
kT = sc.k*300
ni = 1.075*1e+16
N = 1000
Nacc = 1e+23
q = sc.e 
e2 = 11.7*sc.epsilon_0
deltax = ((600*1e-9)/(N-1))
Lx = 1e-7
Ly = 1e-7
mzz = 0.19*sc.m_e

inter = int((N-1)/2)
thermal = kT/q

itr = 20 
centerphi = np.zeros(itr)

def density_phi_cl(gate):
    p_density=np.zeros(N)
    e_density=np.zeros(N)
    p_density[0] = Nacc;p_density[N-1] = ni*ni/Nacc
    e_density[0] = ni*ni/Nacc;e_density[N-1] = Nacc
    density=np.zeros(N)
    phi = np.full(N,0.33374)
    phi[0] = -kT*np.log(Nacc/ni)/q
    phi[N-1] = kT*np.log(Nacc/ni)/q
    #for i in range(itr) :
    GV = (gate)/10.0
    centerphi = 1000
    diff = 1000
    while diff>1e-5 :
        centerphi = phi[int(N/2)]
        e_density = ni*np.exp(q*phi/kT)
        p_density = ni*np.exp(-q*phi/kT)
        p_density[0] = Nacc; p_density[N-1] = ni*ni/Nacc
        e_density[0] = ni*ni/Nacc; e_density[N-1] = Nacc
        phi[0] = -kT*np.log(Nacc/ni)/q
        phi[N-1] = kT*np.log(Nacc/ni)/q

        density = 2*ni*np.sinh(q*phi/kT)

        H = np.zeros((N,N))
        np.fill_diagonal(H,-2*e2)
        np.fill_diagonal(H[1:],e2)
        np.fill_diagonal(H[:,1:],e2)

        H[0][0] = 1
        H[0][1] = 0
        H[N-1][N-1] = 1
        H[N-1][N-2] = 0

        B = np.zeros(N)
        B[:inter] = q*(Nacc+density[:inter])/2.0
        B[inter] = q*density[inter]
        B[inter+1:] = q*(-Nacc+density[inter+1:])/2.0
        B = B*deltax*deltax
        B[0] = phi[0] ; B[N-1] = phi[N-1];

        res = np.zeros(N)
        res = np.matmul(H,phi) - B

        der_B = np.zeros((N,N))
        np.fill_diagonal(der_B[:,:],2*q*(q/kT)*ni*np.cosh(q*phi[:]/kT))
        der_B[inter,inter] = (q/kT)*q*ni*np.cosh(q*phi[inter]/kT)
        np.fill_diagonal(der_B[:inter,:inter],2*q*(q/kT)*ni*np.cosh(q*phi[:inter]/kT))

        der_B = der_B*deltax*deltax

        Jaco = np.zeros((N,N))
        Jaco = H - der_B
        
        phi=phi + slin.solve(Jaco,-res)

        diff = abs(phi[int(N/2)]-centerphi)
    x = np.linspace(0,600*1e-9,num=N,endpoint=True)
    return phi,p_density, e_density

def DD(phi,e_density,p_density):
    VT = q/kT
    Jaco_e = np.zeros((N,N)); Jaco_p = np.zeros((N,N))
    residue_e = np.zeros(N); residue_p = np.zeros(N)

    der_phi = np.zeros(N)
    der_e = np.zeros(N); avg_e = np.zeros(N)
    der_p = np.zeros(N); avg_p = np.zeros(N)

    der_phi[:N-1] = (phi[1:]-phi[:N-1])/deltax
    der_e[:N-1] = (e_density[1:]-e_density[:N-1])/deltax
    der_p[:N-1] = (p_density[1:]-p_density[:N-1])/deltax
    avg_e[:N-1] = (e_density[1:]+e_density[:N-1])/2
    avg_p[:N-1] = (p_density[1:]+p_density[:N-1])/2

    Jn = avg_e*der_phi-VT*der_e
    Jp = -avg_p*der_phi-VT*der_p
    dJn = der_phi/2.0+VT/deltax
    dJn_ = der_phi/2.0-VT/deltax
    dJp = -der_phi/2.0+VT/deltax
    dJp_ = -der_phi/2.0-VT/deltax

    residue_e[1:] = Jn[1:] - Jn[:N-1]
    residue_p[1:] = Jp[1:] - Jp[:N-1]
    np.fill_diagonal(Jaco_e[1:,1:],dJn[1:]-dJn_[:N-1])
    np.fill_diagonal(Jaco_e[:,1:],dJn_[1:])
    np.fill_diagonal(Jaco_e[1:,:],-dJn[1:])
    np.fill_diagonal(Jaco_p[1:,1:],dJp[1:]-dJp_[:N-1])
    np.fill_diagonal(Jaco_p[:,1:],dJp_[1:])
    np.fill_diagonal(Jaco_p[1:,:],-dJp[1:])


    residue_e[0] = e_density[0]-ni*ni/Nacc; residue_p[0] = p_density[0]-Nacc
    residue_e[N-1] = e_density[N-1]-Nacc;  residue_p[N-1] = p_density[N-1]-ni*ni/Nacc
    Jaco_e[0,1] = 0 ; Jaco_e[0,0] = 1
    Jaco_p[0,1] = 0 ; Jaco_p[0,0] = 1
    Jaco_e[N-1,N-2] = 0 ; Jaco_e[N-1,N-1] = 1
    Jaco_p[N-1,N-2] = 0 ; Jaco_p[N-1,N-1] = 1

    dn = slin.solve(Jaco_e,-residue_e)
    dp = slin.solve(Jaco_p,-residue_p)
    e_density = e_density + dn; p_density = p_density + dp
    return phi,e_density, p_density

gate = 0
x = np.linspace(0,600*1e-9,num=N,endpoint=True)
phi, p_density, e_density = density_phi_cl(gate)
old_e = e_density; old_p = p_density
#plt.plot(x, e_density,label=r'$n_{density}$',c='r',marker='o',lw=0,ms=10)
#plt.plot(x, p_density,label=r'$p_{density}$',c='b',marker='o',lw=0,ms=10)

phi, e_density, p_density = DD(phi,e_density,p_density)

plt.plot(x, e_density-old_e,label=r'$n*_{density}-n_{density}$',c='sienna',marker='s',lw=0,ms=10)
plt.plot(x, p_density-old_p,label=r'$p*_{density}-p_{density}$',c='teal',marker='s',lw=0,ms=10)


plt.xlabel("z",fontsize = 18)
plt.ylabel("density difference",fontsize = 18)
#plt.ylabel(r"$density difference$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('symlog')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()


