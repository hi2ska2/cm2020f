import numpy as np
import copy
from scipy.linalg import solve
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import argparse

parser = argparse.ArgumentParser(description='Give Boundary Points')
parser.add_argument('--P0',type=int,help='point of zero value')
parser.add_argument('--P1',type=int,help='point of unit value')
args = parser.parse_args()

vortices = np.loadtxt("HW18_VERTEX.txt")[:,:2]
elements = np.loadtxt("HW18_ELEMENT.txt",dtype=int)-1

def data_processing(Vor,Ele):
    NewVor = np.unique(Vor,axis=0)
    N = NewVor.shape[0]

    NewEle = copy.deepcopy(Ele)
    for i in range(N):
        VorAddr = np.where((Vor==NewVor[i]).all(axis=1))[0]
        for Addr in VorAddr :
            NewEle = np.where(Ele==Addr,i,NewEle)
   
    return NewVor,NewEle

vor,ele = data_processing(vortices,elements)
N = vor.shape[0]

def get_angle(point):
    unit_point = point/np.linalg.norm(point)
    cos = np.dot(unit_point,np.array([1.,0]))

    if point[1] >= 0 : return np.arccos(cos)
    else : return -np.arccos(cos)

def get_length(pos1,pos2):
    return np.sqrt((pos1[0]-pos2[0])**2+(pos1[1]-pos2[1])**2)

def get_circum(tri):
    a,b,c = tri

    d = 2.*(a[0]*(b[1]-c[1])+b[0]*(c[1]-a[1])+c[0]*(a[1]-b[1]))
    ux = (np.dot(a,a)*(b[1]-c[1])+np.dot(b,b)*(c[1]-a[1])+np.dot(c,c)*(a[1]-b[1]))/d
    uy = (np.dot(a,a)*(c[0]-b[0])+np.dot(b,b)*(a[0]-c[0])+np.dot(c,c)*(b[0]-a[0]))/d

    return ux,uy

def get_centroid(tri):
    a,b,c = tri

    ux = (a[0]+b[0]+c[0])/3.
    uy = (a[1]+b[1]+c[1])/3.

    return ux,uy

class Voronoi:
    def __init__(self,ind):
        self.ind = ind
        self.neighbor = set()
        self.boundary = set()

        self.tri_dict = {}
        self.circum_dict = {}
        self.circum_keys = []

    def make_cell(self):
        ele_address = np.where(ele==self.ind)[0]

        for addr in ele_address:
            ele_ = ele[addr]
            tri = np.array([ vor[ele_[j]] for j in range(3) ])

            #tri_plot = mtri.Triangulation(tri[:,0],tri[:,1])
            #plt.triplot(tri_plot,'ro-',ms=12)

            cent = get_centroid(tri)
            circum = get_circum(tri)
            cent_angle = get_angle(cent-vor[self.ind])
            
            #plt.plot(circum[0],circum[1],'g^',ms=12)
            #plt.plot(cent[0],cent[1],'kv',ms=12)
 
            self.tri_dict[cent_angle] = set(ele_)
            self.circum_dict[cent_angle] = circum
            self.neighbor = self.neighbor|set(ele_)

        #plt.plot(vor[self.ind][0],vor[self.ind][1],'bs',ms=12)
        #plt.show()

        self.neighbor.remove(self.ind)
        self.circum_keys = sorted(self.circum_dict)

    def is_boundary(self):
        self.boundary = self.neighbor

        for i in np.arange(len(self.circum_keys))[::-1]:
            inter = self.tri_dict[self.circum_keys[i]]&self.tri_dict[self.circum_keys[i-1]]
            if len(inter) == 2 : 
                self.boundary = self.boundary - inter

        return (self.boundary != set())

def make_H(M1,M0):
    H = np.zeros((N,N))

    for i in range(N):
        cell = Voronoi(i)
        cell.make_cell()
        is_boundary = cell.is_boundary()

        for j in np.arange(len(cell.circum_keys))[::-1]:
            key1 = cell.circum_keys[j]
            key2 = cell.circum_keys[j-1]

            set1 = cell.tri_dict[key1]
            set2 = cell.tri_dict[key2]

            inter = set1&set2
            if len(inter) == 2 :
                inter.remove(i)
                inter = list(inter)
                A = get_length(cell.circum_dict[key1],cell.circum_dict[key2])
                l = get_length(vor[inter[0]],vor[i])
                H[i,inter[0]] += A/l
                H[i,i] -= A/l
            else :
                b1_ = list(cell.boundary&set1)
                if len(b1_) == 1:
                    b1 = b1_[0]
                    b2 = list(cell.boundary&set2)[0]
                else:
                    b1 = b1_[0]
                    b2 = b1_[1]
                A1_proj = np.cross(vor[b1]-vor[i],cell.circum_dict[key1]-vor[i])
                A2_proj = np.cross(vor[b2]-vor[i],cell.circum_dict[key2]-vor[i])
                l1 = get_length(vor[b1],vor[i])
                l2 = get_length(vor[b2],vor[i])
                H[i,b1] += np.abs(A1_proj)/l1**2
                H[i,b2] += np.abs(A2_proj)/l2**2
                H[i,i] -= (np.abs(A1_proj)/l1**2+np.abs(A2_proj)/l2**2)

    M1_arr = np.zeros(N); M1_arr[M1] = 1.
    M0_arr = np.zeros(N); M0_arr[M0] = 1.
    H[M1] = M1_arr; H[M0] = M0_arr

    return H

print("\n Total # of Vertices : "+str(N)+" \n")
M1 = args.P1
M0 = args.P0

H = make_H(M1,M0)
b = np.zeros(N); b[M0] = 0.; b[M1] = 1.

phi = solve(H,b)

triang = mtri.Triangulation(vor[:,0],vor[:,1],ele)
fig = plt.figure()

ax = fig.gca(projection='3d')
cax = ax.plot_trisurf(triang,phi,cmap=cm.jet)
#ax = fig.gca()
#ax.triplot(triang,'ko-',ms=0.5,lw=0.5)

ax.set_xlabel('$x$',fontsize=18)
ax.set_ylabel('$y$',fontsize=18)
ax.set_zlabel('$\phi(x,y)$',fontsize=18)
fig.colorbar(cax)
plt.show()

