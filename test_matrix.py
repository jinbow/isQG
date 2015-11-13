from isqg import *
from numpy import *
from pylab import *
from scipy.io import loadmat, savemat

z = arange(-100,-10,10)
n2 = ones_like(z)*1
f0=1.

nz = len(n2)
    #n2 = r_[1, n2, 1] # Nz->Nz+2 N^2 points including buoyancy at the surface and bottom

zf = r_[z[0]-(z[1]-z[0]), z, z[-1]+(z[-1]-z[-2])] #Nz->Nz+2
    
    #n2 = r_[n2[0], n2]
def mp(a,b):
    return matrix(a)*matrix(b)

zc = r_[zf[0]-0.5*(zf[1]-zf[0]), twopave(zf), zf[-1]+0.5*(zf[-1]-zf[-2])] #Nz+3 \psi points
dzc = diff(zc) # d\psi /dz (Nz+2)
dzf = diff(zf) # (Nz+1)

#nz = n2.size
nz = nz
F = diag(r_[0, f0**2 / n2, 0])
d1 = mp(diag(1. / dzc) , diff(-identity(nz+2),axis=1))
d2 = mp(diag(1./ dzf) , diff(identity(nz+2),axis=0))

eigD2 = mp(d2 , mp(F,d1))

F = diag(r_[f0**2/n2[0], f0**2 / n2, f0**2/n2[-1]])
d2 = mp(diag(r_[1., 1./dzf, 1.]), diff(identity(nz+3),axis=1))
d1 = mp(diag(1. / dzc), diff(-identity(nz+3),axis=0))
sqgD2 = mp(d2,mp(F,d1)); 
sqgD2[-1,-2:]=-1./dzc[-1], 1./dzc[-1]
sqgD2[0,:2]=-1./dzc[0], 1./dzc[0]
