from isqg import *
from numpy import *
from pylab import *
from scipy.io import loadmat, savemat
import popy

N=20.
z = linspace(0,-2000,N)
n2 = exp(z/200)*1e-5
lon = arange(0,7,0.1)
lat = arange(30,37,0.1)
x,y = popy.utils.lonlat2xy(lon,lat)

x = x-x.mean() 
y = y-y.mean()

ssd = 0.1*exp(-(x/5e4)**2-(y/5e4)**2)
ssh = 0.1*exp(-((x-2e4)/7e4)**2-(y/7e4)**2)

useanomaly=True

d = {'z':z,'n2':n2,'lon':lon,'lat':lat,'ssd':ssd,'ssh':ssh,'useanomaly':useanomaly}

savemat('data.mat',d)

contourf(x,y,ssd)
contour(x,y,ssh)
show()