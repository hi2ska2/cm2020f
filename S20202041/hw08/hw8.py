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
#GV=np.linspace(0,1,num=11,endpoint=True)
N = 1000
Nacc = 1e+24
q = sc.e 
e1 = 3.9*sc.epsilon_0
e2 = 11.7*sc.epsilon_0
deltax = ((6.6*1e-9)/N)

inter2 = int((N-1)*(5.8/6.6))
inter1 = int((N-1)*(0.8/6.6))
thermal = kT/q

itr = 20 
centerphi = np.zeros(itr)

for gate in range(11) :
    p_density=np.zeros(N)
    e_density=np.zeros(N)
    density=np.zeros(N)
    phi = np.full(N,0.33374)
    #for i in range(itr) :
    GV = (gate)/10.0
    centerphi = 1000
    diff = 1000
    while diff>1e-5 :
        #GV=0.3
        #centerphi[i] = phi[int(N/2)]
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
        #if(i%3==0) :
    x = np.linspace(0,6.6*1e-9,num=N,endpoint=True)
    plt.plot(x, e_density,label=r'$V_{gate}=$'+str(GV)+'V',c=plt.cm.cool(gate/11),marker='o',lw=0,ms=10)
    plt.xlabel("x",fontsize = 18)
    plt.ylabel(r"$n(r)$",fontsize = 18)
    plt.legend(fontsize=18)
    plt.yscale('symlog')
    plt.tick_params(labelsize=18)
    plt.tight_layout()

plt.show()


'''
plt.plot(range(itr), centerphi,label=r'$V_{gate}=0.3V$',c='r',marker='o',lw=0,ms=10,mfc='none')
plt.xlabel("Iteration",fontsize = 18)
plt.ylabel(r"$\phi(x=3.3nm)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()

plt.plot(x, e_density,label=r'$V_{gate}=0.3V$',c='r',marker='o',lw=0,ms=10,mfc='none')
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$n(r)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
'''

'''
anal_x1 = np.linspace(0,0.8*1e-9,num=1000,endpoint=True)
anal_phi1 = -q*Nacc*5*1e-9*anal_x1/(2*e1)
anal_x2 = np.linspace(0.8*1e-9,5.8*1e-9,num=1000,endpoint=True)
anal_phi2 = (q*Nacc/(2*e2))*(anal_x2-3.3*1e-9)**2-(q*Nacc/2.0)*(4*1e-18/e1+((1e-9*2.5)**2)/e2)
anal_x3 = np.linspace(5.8*1e-9,6.6*1e-9,num=1000,endpoint=True)
anal_phi3 = (q*Nacc*5*1e-9/(2*e1))*(anal_x3-6.6*1e-9)

#plt.plot(x, p_density,c='r',lw=0,marker='o',ms=10,mfc='none')
#plt.plot(x, phi,label=r'$N=$'+str(N),c='r',lw=0,marker='o',ms=10,mfc='none')
#plt.plot(anal_x1, anal_phi1,label='Analytic Solution',c='k',lw=5,marker='o',ms=0,mfc='none')
#plt.plot(anal_x2, anal_phi2,c='k',lw=5,marker='o',ms=0,mfc='none')
#plt.plot(anal_x3, anal_phi3,c='k',lw=5,marker='o',ms=0,mfc='none')
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$n(r)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
'''
