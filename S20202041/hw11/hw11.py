import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

Lz = 5e-9
Ez = sc.hbar*sc.hbar*np.pi*np.pi/(2*0.91*sc.m_e*Lz*Lz)
kT = 300*sc.k
Ef = np.linspace(-0.3*sc.e,0.1*sc.e,10,endpoint=True)

for i in range(np.size(Ef)):
	print("Ef : ",Ef[i])
	err = 1
	density = 0
	n = 1
	while err>1e-5 :
		tmp = density
		density = density + 0.19*sc.m_e*kT*np.log(1+np.exp(-(Ez*n*n-Ef[i])/kT))/(np.pi*sc.hbar*sc.hbar)
		err = abs(tmp-density)
		n = n+1
		print("density : ",density," err : ",err)
	plt.plot(Ef[i], density,label=r'$V_{gate}=1.1$',c='r',marker='o',lw=0,ms=10,mfc='none')

plt.xlabel(r"$E_F(V)$",fontsize = 18)
plt.ylabel(r'$density_{electron}(m^{-2})$',fontsize = 18)
#plt.legend(fontsize=18)
plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()
