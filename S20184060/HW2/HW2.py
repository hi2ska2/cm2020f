import numpy as np
from scipy.linalg import eigh
import scipy.constants as const 
import matplotlib.pyplot as plt

N = 500
L = 5e-9
h = L/(N-1)

H = np.zeros((N-2,N-2))
np.fill_diagonal(H,-2.)
np.fill_diagonal(H[1:],1)
np.fill_diagonal(H[:,1:],1)

E,V = eigh(-H/h**2*const.hbar**2/(2.*0.19*const.m_e))
v0 = np.zeros((1,N-2))
V = np.concatenate((v0,V),axis=0)
V = np.concatenate((V,v0),axis=0)

n = np.arange(1,4)
E_ext = n**2*np.pi**2*const.hbar**2/(2.*0.19*const.m_e*L**2)

for i in range(3):
    print("Numertical "+str(i+1)+"th E : "+str(format(E[i],'8.6e')))
    print("Exact "+str(i+1)+"th E : "+str(format(E_ext[i],'8.6e'))+"\n")

x = np.linspace(0,L,N)
fig = plt.figure(figsize=(10,10))
plt.plot(x,V[:,0]*np.sqrt(N-1),'ro',lw=0.5,ms=5,mfc='None',label='$E_1=$'+str(format(E[0],'8.6e'))+'J')
plt.plot(x,-V[:,1]*np.sqrt(N-1),'bs',lw=0.5,ms=5,mfc='None',label='$E_2=$'+str(format(E[1],'8.6e'))+'J')
plt.plot(x,-V[:,2]*np.sqrt(N-1),'g^',lw=0.5,ms=5,mfc='None',label='$E_3=$'+str(format(E[2],'8.6e'))+'J')

xe = np.linspace(0,L,1000)
plt.plot(xe,np.sin(np.pi*xe/L)*np.sqrt(2),'r-.',lw=0.5,label='$E_1^{ext.}$')
plt.plot(xe,np.sin(2.*np.pi*xe/L)*np.sqrt(2),'b-.',lw=0.5,label='$E_2^{ext.}$')
plt.plot(xe,np.sin(3.*np.pi*xe/L)*np.sqrt(2),'g-.',lw=0.5,label='$E_3^{ext.}$')

plt.xlabel('$x$',fontsize=18)
plt.ylabel(r'$\sqrt{a}\psi_n(x)$',fontsize=18)
legend = plt.legend(loc='upper left', bbox_to_anchor=(1,1),frameon = False, fontsize=12)
legend.set_title('$N=$'+str(N),prop={'size':18})

fig.tight_layout()
plt.show()
