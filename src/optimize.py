"""
Use satellite and ocean in situ data to reconstruct oceanic 3D field.

Created on Mar 20, 2013
@author: Jinbo Wang <jinbow@gmail.com>
"""

import numpy as np
import matplotlib.pyplot as plt
import isqg

# Define longitude and latitude arrays
lon = np.linspace(0,10,50)
lat = np.linspace(30,40,50)

# Convert longitude and latitude arrays to x and y coordinates
x,y = isqg.lonlat2xy(lon,lat)

# Define a Gaussian eddy
radi = 1e5 # a gaussian eddy with 200km radius
gauss = lambda x,y,r: np.exp(- ((x-x.mean())**2 + (y-y.mean())**2) / 2./r**2)
ssd = 0.1 * gauss(x, y, radi)

# Define vertical coordinate array
z = np.linspace(0, -2000, 30)

# Define stratification profile array
n2 = np.ones_like(z) * 1e-5

# Create an isqg_data object with the above parameters
d = isqg.isqg_data(n2=n2, z=z, lon=lon, lat=lat, rho0=1035.)

# Subtract the mean sea surface height (SSH) anomaly from the Gaussian eddy
ssda = ssd - ssd.mean()

# Calculate the mean sea surface height (SSH) anomaly
ssdm = np.ones_like(ssd) * ssd.mean()

# Calculate the mean density field using the stratification profile
rhom = isqg.N2rho(ssdm, isqg.twopave(n2), np.diff(z))

# Set the bottom boundary condition to 'psi=0'
d.bottomboundary = 'psi=0'

# Solve the Quasi-Geostrophic (QG) equations
d.solve_sqg()

# Calculate the SSH from the streamfunction and Rossby deformation radius
d.ssh = np.roll(1.2 * d.psis[0,:,:] * d.f0 / 9.81, 5, axis=1)

# Solve the equations for the streamfunction and velocity potential
d.solve_psii()

# Add the mean SSH anomaly and the mean density field to the QG solution
ssda = ssda + ssdm
d.rhos = d.rhos + rhom

# Plot the results
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))

axs[0, 0].contourf(ssda, levels=np.arange(-0.06,0.11,0.02))
axs[0, 0].set_title('Mean SSH Anomaly')

axs[0, 1].contourf(d.rhos[0,:,:], levels=np.arange(-0.06,0.11,0.02))
axs[0, 1].set_title('Mean Density Field')

axs[1, 0].contourf((ssda-d.rhos[0,:,:])/ssda.max()*100)
axs[1, 0].set_title('Difference between Mean SSH Anomaly and Mean Density Field')

axs[1, 1].contourf(lon, z, d.rhos[:,:,25])
axs[1, 1].set_title('Density Field at Longitude 25')

for ax in axs.flat:
    ax.set(xlabel='Longitude', ylabel='Depth')

plt.tight_layout()
plt.show()

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

axes[0, 0].contourf(d.psit[:, :])
axes[0, 0].set_title('d.psit')

axes[0, 1].contourf(d.psis[0, :, :])
axes[0, 1].set_title('d.psis')

axes[1, 0].contourf(d.psii[0, :, :])
axes[1, 0].set_title('d.psii')

axes[1, 1].contourf(d.psii[0, :, :] + d.psis[0, :, :])
axes[1, 1].set_title('d.psii + d.psis')

plt.tight_layout()
plt.show()
