import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

n_int = 1e+16

N_plus = np.linspace(16,24,num=5000,endpoint=True)
N_plus1 = -np.power(10,N_plus)
N_plus1 = np.flip(N_plus1)
N_plus2 = np.power(10,N_plus)
N_plus = np.append(N_plus1,N_plus2)

T = 300
VT = sc.k*T/sc.e

phi_an = VT*np.arcsinh(0.5*N_plus/n_int)

x = 9
n_plus = np.linspace(16,24,num=x,endpoint=True)
n_plus1 = -np.power(10,n_plus)
n_plus1 = np.flip(n_plus1)
n_plus2 = np.power(10,n_plus)
n_plus = np.append(n_plus1,n_plus2)

#phi_an = VT*np.arcsinh(0.5*n_plus/n_int)

recur =20
president = np.zeros(recur)

phi = np.linspace(-0.7,0.7,num=2*x,endpoint = True)

for j in range(recur) :
    president[j] = np.abs(phi[2*x-1] - phi_an[10000-1])
    if(j%4==0) : plt.plot(n_plus,phi,lw=0,marker='s',label="recur "+str(j),ms=10,zorder=j+1)
    for i in range(2*x) :
        F = n_plus[i]+n_int*np.exp(-phi[i]/VT)-n_int*np.exp(phi[i]/VT)
        partial = -(n_int/VT)*np.exp(-phi[i]/VT)-(n_int/VT)*np.exp(phi[i]/VT)
        delta = -F/partial
        phi[i] = phi[i] + delta





#plt.plot(range(apx_itr+1), x,c='r',lw=0,marker='o',ms=10,mfc='none')
plt.plot(N_plus,phi_an,c='k',label='Analytic',lw=0,marker='o',ms=10,mfc='none',zorder=0)
#plt.plot(n_plus,phi_an,c='r',label='anal',lw=0,marker='o',ms=10,mfc='none',zorder=recur+1)
#plt.plot(anal_x2, anal_phi2,c='k',lw=5,marker='o',ms=0,mfc='none')
#plt.plot(anal_x3, anal_phi3,c='k',lw=5,marker='o',ms=0,mfc='none')
plt.xlabel(r"$N^+$",fontsize = 18)
plt.ylabel(r"$\phi$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.xscale('symlog')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()


plt.plot(range(recur),president,label=r'$N^+=10^{24}m^{-3}$',c='r',lw=0,marker='o',ms=10,mfc='none')
plt.xlabel(r"$iteration$",fontsize = 18)
plt.ylabel(r"$|\phi_{analytic}-\phi_{numerical}|$",fontsize = 18)
plt.legend(fontsize=18)
plt.yscale('log')
#plt.xscale('symlog')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()



