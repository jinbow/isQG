class isqg_data:
    """
    reconstruction ocean 3d field using surface and interior observations
    requirement of initialization:
    N2 on z
    longitude and latitude
    rhos: surface density
    d.rho0: mean density
    d.filterL: the wavelength for the high pass filter
    d.bc1: bottom boundary condition, can be 'psi=0', 'rho=0' (default)
    """

    def __init__(self,n2=0.,z=0.,ssd=0.,ssh=0.,lon=0.,lat=0.,rho0=1035.):
        """set bottom boundary condition
        'rho=0': zero density
        'psi=0': zero pressure
        """
        import sys
        self.bottomboundary = 'rho=0' # zero density at the bottom
        self.filterL = 1e20 # no filtering as default
        self.n2, self.z = n2, z
        self.ssd, self.rho0, self.ssh = ssd, rho0, ssh
        self.lon, self.lat = lon, lat
        self.rstfn = sys.argv[0].replace('py','rst')
        self.filterL = 1e40 # filtered scale
        
        return
    
    def solve(self):
        self.precondition()
        self.solve_sqg()
        self.solve_psii()
        return

    def precondition(self):
        import numpy as np
        print self.rstfn
        from numpy import r_,diff,diag,matrix,identity
        
        
        " --  check n2"
        if self.n2.min()<0:
            print "Error: negative values in N^2"
            return
        if max(self.z) == 0:
            print "error: do not place first N^2 at z=0"
            return
        else:
            self.nz = len(self.n2)
        
        
        " --  set coordinate"
        if np.isscalar(self.lon) and np.isscalar(self.x):
            print "error:  please specify the latitude and longitude"
        elif not np.isscalar(self.lon):
            self.x, self.y = lonlat2xy(self.lon, self.lat)
        else:
            self.x, self.y = np.meshgrid(self.x, self.y)
            
        self.dx = (self.x[:, 1] - self.x[:, 0]).reshape(-1, 1)
        self.dy = (self.y[1, :] - self.y[0, :]).reshape(1, -1)
        self.ny, self.nx = self.x.shape
        " -- initialize parameters"
        self.f0 = lat2f(self.lat.mean())
        
        #set up the vertical difference matrix
        z=self.z;nz=self.nz
        
        zf = r_[z[0]-(z[1]-z[0]), z, z[-1]+(z[-1]-z[-2])] #Nz->Nz+2
            
            #n2 = r_[n2[0], n2]
        def mp(a,b):
            
            return matrix(a)*matrix(b)
        
        zc = r_[zf[0]-0.5*(zf[1]-zf[0]), twopave(zf), zf[-1]+0.5*(zf[-1]-zf[-2])] #Nz+3 \psi points
        dzc = diff(zc) # d\psi /dz (Nz+2)
        dzf = diff(zf) # (Nz+1)
        
        mx = matrix
        f0=self.f0;n2=self.n2
        self.F = mx(diag(r_[f0**2/n2[0], f0**2/n2, f0**2/n2[-1]]))
        self.D1 = mx(mp(diag(1. / dzc) , diff(-identity(nz+2),axis=1)))
        self.D1[[0,-1],:] = 0.
        self.D2 = mx(mp(diag(1./ dzf) , diff(identity(nz+2),axis=0)))
        self.R = mx(np.zeros((nz+2,1)));self.R[-1]=1.
        self.M = self.D2*self.F*self.D1
        self.Rp = self.D2*self.F*self.R
        
        self.zc,self.zf = zc, zf
        self.dzc, self.dzf = dzc, dzf
        
        return 
    
    def solve_sqg(self):
        import scipy.fftpack as fft
        import numpy as np
        from math import pi
        
        self.precondition()
        
        dx = self.dx
        dy = self.dy
        rhos = self.ssd
        
        bhat = fft.fft2( - 9.81 * rhos / self.rho0)  # calculate bouyance 
        ny, nx = rhos.shape
        nz = self.nz
        k = 2 * pi * fft.fftfreq(nx)
        l = 2 * pi * fft.fftfreq(ny)

        ipsihat = np.zeros((nz+3, ny, nx))*complex(0, 0)
        
        Q = np.zeros((nz + 1, 1), dtype='float64'); Q[[0,-1]] = 0.0 # for interior PV, no used in this version
        
        # cutoff value
        ck, cl = 2 * pi / self.filterL, 2 * pi / self.filterL
        # loop through wavenumbers
        bhats = np.zeros_like(bhat)
        for ik in np.arange(k.size):
            for il in np.arange(l.size):
                wv2 = ((k[ik] / dx[il, 0]) ** 2 + (l[il] / dy[0, ik]) ** 2)
                if wv2 > (ck * ck + cl * cl):
                    bhats[il,ik] = bhat[il,ik]
                    right = - bhat[il, ik] / self.f0 * self.Rp
                    left = self.M - wv2 * np.eye(self.nz+1)
                    ipsihat[1:-1, il, ik] = np.linalg.solve(left, right).flatten()
                else:
                    print 'skip k(ik,il)', ik, il, "wv2 = ", wv2
        
        for k in range(1,nz+2):
            ipsihat[k, :, :] = (fft.ifft2(ipsihat[k, :, :]))
        
        if self.bottomboundary == 'psi=0':
            self.psis = np.r_[(np.real(ipsihat)), np.zeros((1,ny,nx))]
        else:
            self.psis = np.real(ipsihat)
            self.psis[0,:,:]= self.psis[1,:,:]
            self.psis[-1,:,:]=self.psis[-2,:,:]-self.dzc[-1]*np.real(fft.ifft2(-bhats))/self.f0
            
        self.rhos = self.psi2rho(self.psis)
        self.us, self.vs = psi2uv(self.lon, self.lat, self.psis)
        
        return
    
    def psi2rho(self,d):
        from numpy import diff
        rho = - self.rho0 * self.f0 / 9.81 * diff(d, axis=0) / \
              self.dzc.reshape(-1, 1, 1)
        return rho
    
    def qgdecomposition(self):
        import numpy as np
        w, vi = np.linalg.eig(self.M)
        vi = vi[:, np.argsort(abs(w))]
        vi = vi / vi[0,:].reshape(1,-1)
        vi = np.r_[vi[0,:],vi,vi[-1,:]]
        w = np.sort(w)
        self.eigenvalue = w
        self.eigenfunction = vi
        return 
    
    def solve_psii(self):
        from numpy import zeros
        self.psit = anomaly(self.lon, self.lat, self.ssh*9.81/self.f0)
        self.ssha = anomaly(self.lon, self.lat, self.ssh)
        
        self.psit = self.ssh*9.81/self.f0
        self.ssha = self.ssh - self.ssh.mean()
        
        psii = self.psit - self.psis[-1,:,:]
        self.qgdecomposition()
        m1, m2 = 0,1
        delta = self.eigenfunction[-1,m1]*self.eigenfunction[0,m2] - self.eigenfunction[0,m1]*self.eigenfunction[-1,m2]
        a =   ( psii * self.eigenfunction[0,m2] + self.psis[0,:,:]*self.eigenfunction[-1,m2]) / delta
        b = - ( psii * self.eigenfunction[0,m1] + self.psis[0,:,:]*self.eigenfunction[-1,m1]) / delta
        self.a = a
        self.b = b
        self.psii = zeros(self.psis.shape)
        for k in range(self.nz+3):
            self.psii[k,:,:] = a * self.eigenfunction[k,m1] + b * self.eigenfunction[k,m2]        
        self.rhoi = self.psi2rho(self.psii)
        self.ui, self.vi = psi2uv(self.lon, self.lat, self.psii)
        
        self.ut, self.vt = self.us+self.ui, self.vs+self.vi
        self.vorticity = vorticity(self.lon,self.lat,self.ut,self.vt)
        
        self.rst = ['.. _'+self.rstfn.replace('.rst','')+':\n\n',
                         'Parameters\n',
                         '==================================\n\n',
                         'bottom boundary condition:'+self.bottomboundary+'\n\n',
                         'size of zc='+str(self.zc.shape)+'\n\n',
                         'size of zf='+str(self.zf.shape)+'\n\n',
                         'size of n2='+str(self.n2.shape)+'\n\n',
                         'size of psis='+str(self.psis.shape)+'\n\n',
                         'size of rhos='+str(self.rhos.shape)+'\n\n',
                         'size of us,vs='+str(self.us.shape)+'\n\n',
                         'size of psii='+str(self.psii.shape)+'\n\n',
                         'size of rhoi='+str(self.rhoi.shape)+'\n\n',
                         'size of ui,vi='+str(self.ui.shape)+'\n\n',
                         '================  ================  ================  ================\n',
                         ' i                  zc                     zf             :math:`N_2` \n',
                         '================  ================  ================  ================\n',]
        for i in range(min(self.zc.size,self.zf.size,self.n2.size)):
            self.rst.append('%16i%16.2f%16.2f%16.2e\n'%(i,self.zc[i],self.zf[i],self.n2[i]))
                    
        self.rst.append('%16i%16.2f%16.2f%16.2e\n'%(i+1,self.zc[i+1],0,0))

        self.rst.append('================  ================  ================  ================\n\n',)
        return
    
    def plotgrid(self):
        """
        show the grid structure on a figure
        """
        from pylab import subplot
        from numpy import zeros
        ax = subplot(111)
        ax.plot(zeros(4), self.zf[:4],'ko',markerfacecolor='w',markersize=10)
        for i in range(4):
            ax.text(0.4, self.zf[i]*1.03, "$zf_%1i$"%i,fontsize=20) 
            ax.text(-0.9, self.zc[i]*1.02, "$zc_%1i$"%i,fontsize=20)
        ax.plot(zeros(4), self.zc[:4], 'k*')
        
        ax.plot(zeros(4),self.zf[-4:]-self.zf[-8],'ko',markerfacecolor='w',markersize=10)
        ax.plot(zeros(4), self.zc[-4:]-self.zf[-8], 'k*')
        
        ax.axis([-5,5,self.zf[-1]-self.zf[-8]-50,50])
        ax.set_yticks([])
        ax.set_xticks([])
        return
    


def fit2Dsurf(x,y,p):
    """
      given y0=f(t0), find the best fit
      p = a + bx + cy + dx**2 + ey**2 + fxy
      and return a,b,c,d,e,f
    """
    from scipy.optimize import leastsq
    def fit2poly(c,x,y,p):
        a,b,c,d,e,f=c
        err = p - (a + b*x + c*y + d*x**2 + e*y**2 + f*x*y)
        return err
    x=x.flatten()
    y=y.flatten()
    p=p.flatten()
    c = [p.mean(),1e-3,1e-3,1e-6,1e-6,1e-6]
    coef = leastsq(fit2poly,c,args=(x,y,p))
    return coef

def fit2exp(x,y,method='exp'):
    """
      given y0=f(t0), find the best fit
      p = a + b*exp(c*z)
      and return a,b,c,d,e,f
    """
    from scipy.optimize import leastsq
    import numpy as np

    def fit2poly(cc,x,y):
        if method=='exp':
            err = y - ( cc[0]*np.exp(cc[1]*x))
        elif method=='tanh':
            err = y - ( cc[0]*(1. - np.tanh(cc[1]*x)**2))
        return err

    x=x.flatten()
    y=y.flatten()
    c = [1e-5,1./200.]
    coef = leastsq(fit2poly,c,args=(x,y))
    return coef


def anomaly(lon,lat,var):
    from numpy import meshgrid,zeros,arange
    if lon.ndim == 1:
        x1,y1=meshgrid(lon,lat)
    else:
        x1,y1 = lon, lat
    if var.ndim == 2:
        coef = fit2Dsurf(x1,y1,var)[0]
        varm = fit2poly(coef,x1,y1)
        return var-varm
    elif var.ndim == 3:
        tmp = zeros(var.shape)
        for i in arange(var.shape[0]):
            coef = fit2Dsurf(x1,y1,var[i,:,:])[0]
            varm = fit2poly(coef,x1,y1)
            tmp[i,:,:] = var[i,:,:]-varm
        return tmp

def fit2poly(c,x,y):
    a,b,c,d,e,f=c
    fit = (a + b*x + c*y + d*x**2 + e*y**2 + f*x*y)
    return fit

def lat2f(d):
    """ Calculate Coriolis parameter from latitude
    d: latitudes, 1- or 2-D
    """
    from math import pi, sin
    return 2.0*0.729e-4*sin(d*pi/180.0)

def lonlat2xy(lon,lat):
    from pylab import meshgrid,cos,pi
    r = 6371.e3
    lon = lon-lon[0]
    if lon.ndim == 1:
        lon,lat = meshgrid(lon,lat)
    x = 2*pi*r*cos(lat*pi/180.)*lon/360.
    y = 2*pi*r*lat/360.
    return x,y

def vorticity(lon,lat,u,v):
    tmp, dudy = gradxy(lon,lat,u)
    dvdx, tmp = gradxy(lon,lat,v)
    del tmp
    return dvdx-dudy

def gradxy(lon,lat,d):
    """ calculate the x and y gradients of d
    d can be of two or three dimensions.
    the calculated velocity fields share the same dimensions with d,
    i.e. the velocity fields are shifted.
    """
    from numpy import c_,r_,diff,zeros
    x,y = lonlat2xy(lon,lat)
    def uv(x,y,d):
        x = c_[x, x[:,-2]]
        y = r_[y, y[-2,:].reshape(1,-1)]
        
        sshx = c_[d, d[:,-1]]
        v = diff(sshx,axis=1)/diff(x,axis=1)
        v[:,-1] = v[:,-2]

        sshy = r_[d, d[-1,:].reshape(1,-1)]
        u = diff(sshy,axis=0)/diff(y,axis=0)
        u[-2,:]=u[-1,:]
        return u,v
    if d.ndim == 2:
        return uv(x,y,d)
    elif d.ndim == 3:
        u, v = zeros(d.shape), zeros(d.shape)
        for i in range(d.shape[0]):
            ut,vt = uv(x,y,d[i,:,:])
            v[i,:,:] = vt
            u[i,:,:] = ut
        return u, v
    
def psi2uv(lon,lat,d):
    """calculate horizontal velocity from streamfunction
    """
    u, v = gradxy(lon, lat, d)
    u = -1 * u
    return u, v

def ssh2uv(lon,lat,ssh):
    """ calculate geostrophic velocity from SSH data
    lon[nx], 
    lat[ny]
    ssh[ny,nx] """
    from pylab import meshgrid
    if len(lon.shape) == 1:
        lon,lat = meshgrid(lon,lat)
    psi = 9.81 / lat2f(lat) * ssh
    u, v = gradxy(lon,lat,psi)
    u = -1 * u
    return u, v

def twopave(x):
    return (x[0:-1]+x[1:])/2.

def N2rho(rhos,n2,dz,rho0=1035., g = 9.81, direction='top2bottom'):
    """
    Calculate the mean density field from a surface density field and a stratification profile,
    :math:`\rho = - {\rho_0 over g} \int_z N^2 dz`
    """
    import numpy as np
    nz = n2.size; ny,nx = rhos.shape
    rho = np.zeros((nz+1,ny,nx))
    rho[0,:,:] = rhos
    dz = abs(dz)
    for i in range(1,nz+1):
        if direction=='top2bottom':
            rho[i,:,:] = rho[i-1,:,:] + dz[i-1]*n2[i-1]*rho0/g
        else:
            rho[i,:,:] = rho[i-1,:,:] - dz[i-1]*n2[i-1]*rho0/g
    return rho