import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as sc
import scipy.linalg as slin
import scipy.integrate as inte
import fractions

#np.set_printoptions(formatter={'all':lambda x: str(fractions.Fraction(x).limit_denominator())})


N = 100

e1 = 3.9
e2 = 22
H = np.zeros((N,N))

np.fill_diagonal(H,-2*e2)
np.fill_diagonal(H[1:],e2)
np.fill_diagonal(H[:,1:],e2)

inter = int((N-1)*(0.5/2.4))
np.fill_diagonal(H[0:inter+1],-2*e1)
np.fill_diagonal(H[1:inter+1],e1)
np.fill_diagonal(H[:,1:inter+1],e1)


H[inter][inter] = -e2-e1


H[0][0] = 1
H[0][1] = 0
H[N-1][N-1] = 1
H[N-1][N-2] = 0


B = np.zeros(N)
B[0] = 0
B[N-1] = 1

X = np.linalg.solve(H,B)
x = np.linspace(0,2.4*1e-9,num=N,endpoint=True)


plt.plot(x, X,label=r'$N=$'+str(N),c='r',lw=0,marker='o',ms=10,mfc='none')
plt.xlabel("x",fontsize = 18)
plt.ylabel(r"$\phi(x)$",fontsize = 18)
plt.legend(fontsize=18)
#plt.yscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()

dy_1 = X[1]
dx = x[1]
dy_2 = X[N-1] - X[N-2]

slope1 = dy_1/dx
slope2 = dy_2/dx

e0 = sc.epsilon_0
C = 1e+9*e0*e1*2200/1841

print("Analytic C = ",C)


numb= np.zeros(4)
C_nu = np.zeros(4)
for i in range(4) :
    N = 10**(i+1)
    numb[i]=N
    e1 = 3.9
    e2 = 22
    H = np.zeros((N,N))
    np.fill_diagonal(H,-2*e2)
    np.fill_diagonal(H[1:],e2)
    np.fill_diagonal(H[:,1:],e2)
    inter = int((N-1)*(0.5/2.4))
    np.fill_diagonal(H[0:inter+1],-2*e1)
    np.fill_diagonal(H[1:inter+1],e1)
    np.fill_diagonal(H[:,1:inter+1],e1)
    H[inter][inter] = -e2-e1
    H[0][0] = 1
    H[0][1] = 0
    H[N-1][N-1] = 1
    H[N-1][N-2] = 0
    B = np.zeros(N)
    B[0] = 0
    B[N-1] = 1
    X = np.linalg.solve(H,B)
    x = np.linspace(0,2.4*1e-9,num=N,endpoint=True)
    dy_1 = X[1]
    dx = x[1]
    slope1 = dy_1/dx
    e0 = sc.epsilon_0
    C_nu[i]=slope1*e1*e0
    print("N = ",N)
    print("C_numerical = ", C_nu[i])


plt.plot(numb, abs(C_nu-C),label=r'$N=$'+str(N),c='r',lw=0,marker='o',ms=10,mfc='none')
plt.xlabel("N",fontsize = 18)
plt.ylabel(r"$|C-C_{numerical}|$",fontsize = 18)
#plt.legend(fontsize=18)
plt.yscale('log')
plt.xscale('log')
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.show()


