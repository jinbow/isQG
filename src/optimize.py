'''
Use satellite and ocean in situ data to reconstruct oceanic 3D field.

Created on Mar 20, 2013

@author: Jinbo Wang <jinbow@gmail.com>
'''


if __name__ == '__main__':
    import numpy as np
    import pylab as plt
    import isqg
    
    lon = np.linspace(0,10,50)
    lat = np.linspace(30,40,50)
    x,y = isqg.lonlat2xy(lon,lat)
    radi = 1e5 # a gaussian eddy with 200km radius
    gauss = lambda x,y,r: np.exp(- ( 
                     ( x-x.mean() )**2 + ( y-y.mean() )**2 
                       )  /2./r**2 )
    ssd = 0.1*gauss(x,y,radi)
    z = np.linspace(0,-2000,30)
    n2 = np.ones_like(z)* 1e-5
    
    d = isqg.isqg_data(n2=n2,z=z,lon=lon,lat=lat,rho0=1035.)
    ssda=ssd-ssd.mean()
    ssdm = np.ones_like(ssd)*ssd.mean()
    rhom = isqg.N2rho(ssdm, isqg.twopave(n2), np.diff(z))
    d.bottomboundary='psi=0'
    d.solve_sqg()
    d.ssh = np.roll(1.2*d.psis[0,:,:]*d.f0/9.81,5,axis=1)
    d.solve_psii()
            
    ssda = ssda+ssdm
    d.rhos = d.rhos+rhom
    
    plt.figure()
    plt.subplot(221)
    plt.contourf(ssda,levels=np.arange(-0.06,0.11,0.02))
    plt.colorbar()
    plt.subplot(222)
    plt.contourf(d.rhos[0,:,:],levels=np.arange(-0.06,0.11,0.02))
    plt.colorbar()
    plt.subplot(223)
    plt.contourf((ssda-d.rhos[0,:,:])/ssda.max()*100)
    plt.colorbar()
    plt.subplot(224)
    plt.contourf(lon,z,d.rhos[:,:,25])
    plt.colorbar()
    plt.figure()
    plt.subplot(221)
    plt.contourf(d.psit[:,:])
    plt.colorbar()
    plt.contour(d.psit)
    plt.subplot(222)
    plt.contourf(d.psis[0,:,:])
    plt.colorbar()
    plt.contour(d.psit)
    plt.subplot(223)
    plt.contourf(d.psii[0,:,:])
    plt.colorbar()
    plt.contour(d.psit)
    plt.subplot(224)
    plt.contourf(d.psii[0,:,:]+d.psis[0,:,:])
    plt.colorbar()
    plt.contour(d.psit)
    plt.show()
    
