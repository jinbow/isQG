from isqg import *
from numpy import *
from pylab import *
from scipy.io import loadmat, savemat

d = isqg_data()
N=50.
d.z = linspace(-500,-10,N)
d.n2 = exp(d.z/200)*1e-5
d.n2 = ones_like(d.z)*1e-5
d.lon = arange(0,7,0.5)
d.lat = arange(30,37,0.5)
d.precondition()
d.ssd = 0.01*cos(2*pi*d.x/d.x.max())*cos(2*pi*d.x.T/d.x.max())
x,y=d.x,d.y
d.ssd = 0.1*exp(-((x-x.mean())/9e4)**2-((y-y.mean())/9e4)**2)
d.rho0 = 1035.
d.filterL = 1e8
d.bottomboundary = 'rho=0'

d.solve_sqg()
 
figure(0)
subplot(221)
contourf(d.ssd,10);colorbar()
subplot(222)
contourf(d.rhos[-1,:,:],10);colorbar()
subplot(223)
ax=pcolor(d.rhos[-1,:,:]-d.ssd);cb=colorbar(ax);
subplot(224)
contourf(d.lon, d.zf, d.rhos[:,d.ny/2,:],10);colorbar()
 
figure(1)
 
d.plotgrid()
 
show()
