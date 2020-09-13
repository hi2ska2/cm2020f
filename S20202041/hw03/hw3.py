import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte

N = 60

H = np.zeros((N,N))
np.fill_diagonal(H,-2)
np.fill_diagonal(H[1:],1)
np.fill_diagonal(H[:,1:],1)

H = H*N*N

H[0][0] = 1
H[0][1] = 0
H[N-1][N-1] = 1
H[N-1][N-2] = 0

B = np.zeros(N)
B[0] = 1
B[N-1] = -1

X = np.linalg.solve(H,B)
x = np.linspace(0,1,num=N,endpoint=True)

plt.plot(x,-2*x+1,label=r'$y=-2x+1$',c='k',lw=5,ms=0)
plt.plot(x, X,label=r'$N=6$',c='r',lw=0,marker='o',ms=10,mfc='none')
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$\phi(x)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()


#Problem 2, delta function

delta = np.zeros(N)
delta[(int)(N/2)-1] = N-1
delta[(int)(N/2)] = N-1

#Plotting Delta function, square function.
for i in range(3):
    if i == 0 : 
        N_delta = 6
    elif i == 1 :
        N_delta = 18
    else :
        N_delta = 30
    B1 = np.zeros(N_delta)
    B1[(int)((N_delta)/2)-1] = N_delta-1
    B1[(int)((N_delta)/2)] = N_delta-1
    x1 = np.linspace(0,1,num=N_delta)
    plt.plot(x1, B1,label='N = '+str(N_delta),marker='o',lw=0,ms=10,mfc='none')

plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$\delta(x)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()


X = np.linalg.solve(H,delta)
N1=100
x2 = np.linspace(0,1,num=N1)

plt.plot(x2[:int(N1/2)],-x2[:int(N1/2)]/2,c='k',lw=5,marker='o',ms=0,mfc='none')
plt.plot(x2[int(N1/2):],x2[int(N1/2):]/2-1/2,label='Analytic sol',c='k',lw=5,marker='o',ms=0,mfc='none')
plt.plot(x, X,label='N = '+str(N),c='r',lw=0,marker='o',ms=10,mfc='none')
plt.xlim([0,1])
plt.ylim([-0.25,0])
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$\phi(x)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()


