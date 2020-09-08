import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte

E_gap = np.zeros(3)
n = 1
L = 5*1e-9
m = 0.91*sc.m_e
hbar = sc.hbar

for i in range(3) :
        N = 5*10**i
        x = np.arange(N)
        x = L*x/(N-1)
        dx = x[1]
        new_x = np.delete(x,[0,N-1])

        #numerical solution
        H = np.zeros((N-2,N-2))
        np.fill_diagonal(H,-2)
        np.fill_diagonal(H[1:],1)
        np.fill_diagonal(H[:,1:],1)

        E, psi = slin.eigh(-hbar*hbar*H/(2*m*dx*dx))
        norm = inte.simps(psi[:,n]*psi[:,n],new_x)
        if(i==2) :
                psi[:,n] = psi[:,n]/np.sqrt(norm)
        else :
                psi[:,n] = -psi[:,n]/np.sqrt(norm)
        plt.plot(new_x,psi[:,n],label='numerical',c='r',lw=0,marker='o',ms=10,mfc='none')
        print("N = ",N)
        print("Numerical energy = ", E[n])
        #analytic solution
        x = np.arange(1000)
        x = L*x/999
        psi_a = np.sqrt(2/L)*np.sin((n+1)*np.pi*x/L)
        plt.plot(x,psi_a,label='analytic',c='k',lw=5,marker='d',ms=0,mfc='none')

        E_anal = (n+1)*(n+1)*np.pi*np.pi*hbar*hbar/(2*m*L*L)
        print("Analytic energy = ", E_anal)
        plt.xlabel("x",fontsize = 18)
        plt.ylabel("y",fontsize = 18)
        plt.legend(fontsize=18)
        plt.tick_params(labelsize=18)
        plt.tight_layout()
        #plt.savefig('hw2.pdf')
        plt.show()

        E_gap[i]=abs(E_anal- E[n])


plt.plot([5,50,500], E_gap,label=r'$E_{numerical} - E_{analytic}$',c='r',lw=0,marker='o',ms=10,mfc='none')
plt.xlabel("N",fontsize = 18)
plt.ylabel(r"$|E_{num}-E_{ana}|$",fontsize = 18)
#plt.legend(fontsize=18)
plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
                              
