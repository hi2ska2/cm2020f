import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions

#np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})


N = 1000
Nacc = 10**25
q = sc.e 

e1 = 3.9*sc.epsilon_0
e2 = 11.7*sc.epsilon_0

H = np.zeros((N,N))

np.fill_diagonal(H,-2*e1)
np.fill_diagonal(H[1:],e1)
np.fill_diagonal(H[:,1:],e1)

inter2 = int((N-1)*(5.8/6.6))
np.fill_diagonal(H[0:inter2+1],-2*e2)
np.fill_diagonal(H[1:inter2+1],e2)
np.fill_diagonal(H[:,1:inter2+1],e2)
H[inter2][inter2] = -e2-e1

inter1 = int((N-1)*(0.8/6.6))
np.fill_diagonal(H[0:inter1+1],-2*e1)
np.fill_diagonal(H[1:inter1+1],e1)
np.fill_diagonal(H[:,1:inter1+1],e1)
H[inter1][inter1] = -e2-e1

H[0][0] = 1
H[0][1] = 0
H[N-1][N-1] = 1
H[N-1][N-2] = 0


print(H)

print("inter1 : ", inter1)
print("inter2 : ", inter2)

B = np.zeros(N)
B[inter1] = q*Nacc/2.0
B[inter1+1:inter2] = q*Nacc
B[inter2] = q*Nacc/2.0
B = B*((6.6*1e-9)/N)**2

print(B)


X = np.linalg.solve(H,B)
x = np.linspace(0,6.6*1e-9,num=N,endpoint=True)

anal_x1 = np.linspace(0,0.8*1e-9,num=1000,endpoint=True)
anal_phi1 = -q*Nacc*5*1e-9*anal_x1/(2*e1)
anal_x2 = np.linspace(0.8*1e-9,5.8*1e-9,num=1000,endpoint=True)
anal_phi2 = (q*Nacc/(2*e2))*(anal_x2-3.3*1e-9)**2-(q*Nacc/2.0)*(4*1e-18/e1+((1e-9*2.5)**2)/e2)
anal_x3 = np.linspace(5.8*1e-9,6.6*1e-9,num=1000,endpoint=True)
anal_phi3 = (q*Nacc*5*1e-9/(2*e1))*(anal_x3-6.6*1e-9)

plt.plot(x, X,label=r'$N=$'+str(N),c='r',lw=0,marker='o',ms=10,mfc='none')
plt.plot(anal_x1, anal_phi1,label='Analytic Solution',c='k',lw=5,marker='o',ms=0,mfc='none')
plt.plot(anal_x2, anal_phi2,c='k',lw=5,marker='o',ms=0,mfc='none')
plt.plot(anal_x3, anal_phi3,c='k',lw=5,marker='o',ms=0,mfc='none')
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$\phi(x)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()

