import sys
import logging
from typing import Union, Tuple, List, Optional
import numpy as np
import scipy.fftpack as fft
import matplotlib.pyplot as plt
from numpy import ndarray
from scipy import sparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)



class isqg_data:
    """
    Class for reconstructing a 3D ocean field using surface and interior observations.

    This class initializes the necessary parameters for solving the Surface Quasi-Geostrophic (SQG)
    and interior streamfunction equations, given stratification (N^2), coordinates, and boundary conditions.

    Attributes:
        bottomboundary (str): Bottom boundary condition ('rho=0' for zero density, 'psi=0' for zero pressure).
        filterL (float): Filter wavelength for wavenumber cutoff (m).
        n2 (ndarray): Stratification profile (N^2) on z-levels (s^-2).
        z (ndarray): Vertical coordinates (m).
        ssd (ndarray): Surface density anomaly (kg/m^3).
        rho0 (float): Mean density (kg/m^3).
        ssh (ndarray): Sea surface height anomaly (m).
        lon (ndarray): Longitude coordinates (degrees).
        lat (ndarray): Latitude coordinates (degrees).
        rstfn (str): Filename for RST documentation output.
        x (ndarray): X-coordinates in Cartesian plane (m).
        y (ndarray): Y-coordinates in Cartesian plane (m).
        dx (ndarray): X-grid spacing (m).
        dy (ndarray): Y-grid spacing (m).
        nx (int): Number of x-grid points.
        ny (int): Number of y-grid points.
        nz (int): Number of z-grid points.
        f0 (float): Coriolis parameter (s^-1).
        zc (ndarray): Vertical coordinates for cell centers (m).
        zf (ndarray): Vertical coordinates for cell faces (m).
        dzc (ndarray): Vertical spacing for cell centers (m).
        dzf (ndarray): Vertical spacing for cell faces (m).
        F (ndarray): Vertical stratification matrix.
        D1 (ndarray): First-order vertical difference matrix.
        D2 (ndarray): Second-order vertical difference matrix.
        R (ndarray): Boundary condition vector.
        M (ndarray): Vertical operator matrix.
        Rp (ndarray): Boundary condition operator.
        psis (ndarray): Surface streamfunction (m^2/s).
        rhos (ndarray): Density anomaly field (kg/m^3).
        us (ndarray): Surface zonal velocity (m/s).
        vs (ndarray): Surface meridional velocity (m/s).
        eigenvalue (ndarray): Eigenvalues from QG decomposition.
        eigenfunction (ndarray): Eigenfunctions from QG decomposition.
        psit (ndarray): Total streamfunction (m^2/s).
        ssha (ndarray): Sea surface height anomaly (m).
        psii (ndarray): Interior streamfunction (m^2/s).
        rhoi (ndarray): Interior density anomaly (kg/m^3).
        ui (ndarray): Interior zonal velocity (m/s).
        vi (ndarray): Interior meridional velocity (m/s).
        ut (ndarray): Total zonal velocity (m/s).
        vt (ndarray): Total meridional velocity (m/s).
        vorticity (ndarray): Vorticity field (s^-1).
        rst (List[str]): RST documentation strings.
    """

    def __init__(
        self,
        n2: Union[float, ndarray] = 0.0,
        z: Union[float, ndarray] = 0.0,
        ssd: Union[float, ndarray] = 0.0,
        ssh: Union[float, ndarray] = 0.0,
        lon: Union[float, ndarray] = 0.0,
        lat: Union[float, ndarray] = 0.0,
        rho0: float = 1035.0
    ) -> None:
        """
        Initialize the ISQGData object with default or provided parameters.

        Args:
            n2 (Union[float, ndarray]): Stratification profile (N^2) on z-levels (s^-2).
            z (Union[float, ndarray]): Vertical coordinates (m).
            ssd (Union[float, ndarray]): Surface density anomaly (kg/m^3).
            ssh (Union[float, ndarray]): Sea surface height anomaly (m).
            lon (Union[float, ndarray]): Longitude coordinates (degrees).
            lat (Union[float, ndarray]): Latitude coordinates (degrees).
            rho0 (float): Mean density (kg/m^3). Defaults to 1035.0.

        Raises:
            ValueError: If input parameters are invalid.
        """
        self.bottomboundary: str = 'rho=0'  # Zero density at bottom
        self.filterL: float = 1e20  # Default filter wavelength (no filtering)
        self.n2: Union[float, ndarray] = np.asarray(n2) if not isinstance(n2, float) else n2
        self.z: Union[float, ndarray] = np.asarray(z) if not isinstance(z, float) else z
        self.ssd: Union[float, ndarray] = np.asarray(ssd) if not isinstance(ssd, float) else ssd
        self.rho0: float = rho0
        self.ssh: Union[float, ndarray] = np.asarray(ssh) if not isinstance(ssh, float) else ssh
        self.lon: Union[float, ndarray] = np.asarray(lon) if not isinstance(lon, float) else lon
        self.lat: Union[float, ndarray] = np.asarray(lat) if not isinstance(lat, float) else lat
        self.rstfn: str = sys.argv[0].replace('.py', '.rst')
        logger.debug(f"Initialized ISQGData with rstfn: {self.rstfn}")

    def solve(self) -> None:
        """
        Precondition the data and solve the SQG and PSI equations.

        Raises:
            RuntimeError: If preconditioning or solving fails.
        """
        try:
            self.precondition()
            self.solve_sqg()
            self.solve_psii()
        except Exception as e:
            logger.error(f"Failed to solve equations: {e}")
            raise RuntimeError(f"Solve error: {e}")

    def precondition(self) -> None:
        """
        Set up the vertical difference matrix and initialize parameters.

        Validates the stratification profile (N^2) and coordinates, and constructs
        the vertical difference matrices for the SQG solver.

        Raises:
            ValueError: If N^2 is negative or coordinates are invalid.
            RuntimeError: If matrix construction fails.
        """
        try:
            # Validate N^2
            if not np.isscalar(self.n2) and np.any(self.n2 < 0):
                raise ValueError("Negative values detected in N^2")
            if not np.isscalar(self.z) and np.max(self.z) == 0:
                raise ValueError("First N^2 cannot be at z=0")
            
            self.nz: int = len(self.n2) if not np.isscalar(self.n2) else 1
            logger.debug(f"Number of vertical levels (nz): {self.nz}")

            # Set coordinates
            if np.isscalar(self.lon) or np.isscalar(self.lat):
                raise ValueError("Latitude and longitude must be arrays")
            self.x, self.y = lonlat2xy(self.lon, self.lat)
            self.dx = (self.x[:, 1] - self.x[:, 0]).reshape(-1, 1)
            self.dy = (self.y[1, :] - self.y[0, :]).reshape(1, -1)
            self.ny, self.nx = self.x.shape
            logger.debug(f"Grid dimensions: nx={self.nx}, ny={self.ny}")

            # Initialize Coriolis parameter
            self.f0 = lat2f(np.mean(self.lat))

            # Set up vertical difference matrix
            z = np.asarray(self.z) if not np.isscalar(self.z) else np.array([self.z])
            zf = np.r_[z[0] - (z[1] - z[0]), z, z[-1] + (z[-1] - z[-2])]  # Nz+2
            zc = np.r_[
                zf[0] - 0.5 * (zf[1] - zf[0]),
                twopave(zf),
                zf[-1] + 0.5 * (zf[-1] - zf[-2])
            ]  # Nz+3
            dzc = np.diff(zc)  # Nz+2
            dzf = np.diff(zf)  # Nz+1

            # Matrix operations (use sparse matrices for large nz)
            F = sparse.diags(np.r_[self.f0 ** 2 / self.n2[0], self.f0 ** 2 / self.n2, self.f0 ** 2 / self.n2[-1]])
            D1 = sparse.diags(1. / dzc) @ np.diff(-np.identity(self.nz + 2), axis=1)
            D1[[0, -1], :] = 0.0
            D2 = sparse.diags(1. / dzf) @ np.diff(np.identity(self.nz + 2), axis=0)
            R = np.zeros((self.nz + 2, 1))
            R[-1] = 1.0
            M = D2 @ F @ D1
            Rp = D2 @ F @ R

            # Store attributes
            self.zc, self.zf = zc, zf
            self.dzc, self.dzf = dzc, dzf
            self.F, self.D1, self.D2 = F, D1, D2
            self.R, self.M, self.Rp = R, M, Rp
            logger.info("Preconditioning completed")
        except Exception as e:
            logger.error(f"Precondition failed: {e}")
            raise RuntimeError(f"Precondition error: {e}")

    def solve_sqg(self) -> None:
        """
        Solve the Surface Quasi-Geostrophic (SQG) equation.

        Computes the streamfunction and derived fields (density, velocities) from
        the surface density anomaly.

        Raises:
            RuntimeError: If the SQG solver fails.
        """
        from math import pi
        try:
            self.precondition()
            dx, dy = self.dx, self.dy
            rhos = np.asarray(self.ssd)
            ny, nx = rhos.shape
            nz = self.nz

            # Compute buoyancy in Fourier space
            bhat = fft.fft2(-9.81 * rhos / self.rho0)
            k = 2 * pi * fft.fftfreq(nx)
            l = 2 * pi * fft.fftfreq(ny)
            ipsihat = np.zeros((nz + 3, ny, nx), dtype=complex)

            # Initialize interior PV (not used in this version)
            Q = np.zeros((nz + 1, 1), dtype='float64')
            Q[[0, -1]] = 0.0

            # Wavenumber cutoff
            ck, cl = 2 * pi / self.filterL, 2 * pi / self.filterL
            bhats = np.zeros_like(bhat)
            
            # Solve for each wavenumber
            for ik in range(k.size):
                for il in range(l.size):
                    wv2 = ((k[ik] / dx[il, 0]) ** 2 + (l[il] / dy[0, ik]) ** 2)
                    if wv2 > (ck * ck + cl * cl):
                        bhats[il, ik] = bhat[il, ik]
                        right = -bhat[il, ik] / self.f0 * self.Rp
                        left = self.M - wv2 * np.eye(self.nz + 1)
                        ipsihat[1:-1, il, ik] = np.linalg.solve(left, right).flatten()
                    else:
                        logger.debug(f"Skipping k({ik},{il}), wv2={wv2}")

            # Inverse FFT
            for k in range(1, nz + 2):
                ipsihat[k, :, :] = fft.ifft2(ipsihat[k, :, :])

            # Apply boundary conditions
            if self.bottomboundary == 'psi=0':
                self.psis = np.r_[np.real(ipsihat), np.zeros((1, ny, nx))]
            else:
                self.psis = np.real(ipsihat)
                self.psis[0, :, :] = self.psis[1, :, :]
                self.psis[-1, :, :] = (
                    self.psis[-2, :, :] - self.dzc[-1] * np.real(fft.ifft2(-bhats)) / self.f0
                )

            # Compute derived fields
            self.rhos = self.psi2rho(self.psis)
            self.us, self.vs = psi2uv(self.lon, self.lat, self.psis)
            logger.info("SQG solution completed")
        except Exception as e:
            logger.error(f"SQG solver failed: {e}")
            raise RuntimeError(f"SQG solver error: {e}")

    def psi2rho(self, d: ndarray) -> ndarray:
        """
        Convert streamfunction anomaly to density anomaly.

        Args:
            d (ndarray): Streamfunction anomaly (m^2/s).

        Returns:
            ndarray: Density anomaly (kg/m^3).

        Raises:
            ValueError: If input shape is invalid.
        """
        try:
            if d.ndim != 3:
                raise ValueError(f"Expected 3D streamfunction, got shape {d.shape}")
            rho = -self.rho0 * self.f0 / 9.81 * np.diff(d, axis=0) / self.dzc.reshape(-1, 1, 1)
            return rho
        except Exception as e:
            logger.error(f"Failed to convert psi to rho: {e}")
            raise ValueError(f"psi2rho error: {e}")

    def qgdecomposition(self) -> None:
        """
        Perform Quasi-Geostrophic (QG) decomposition.

        Computes eigenvalues and eigenfunctions of the vertical operator matrix.

        Raises:
            RuntimeError: If eigenvalue decomposition fails.
        """
        try:
            # Convert sparse matrix to dense if necessary
            M_dense = self.M.toarray() if sparse.issparse(self.M) else self.M
            w, vi = np.linalg.eig(M_dense)
            idx = np.argsort(np.abs(w))
            w = w[idx]
            vi = vi[:, idx]
            
            # Normalize eigenfunctions
            vi = vi / vi[0, :].reshape(1, -1)
            # Extend eigenfunctions to include boundary points
            vi = np.r_[vi[0:1, :], vi, vi[-1:, :]]
            
            self.eigenvalue = w
            self.eigenfunction = vi
            logger.info("QG decomposition completed")
        except Exception as e:
            logger.error(f"QG decomposition failed: {e}")
            raise RuntimeError(f"QG decomposition error: {e}")

    def solve_psii(self) -> None:
        """
        Solve for the interior streamfunction anomaly.

        Computes the interior streamfunction, density, velocities, and vorticity
        by decomposing the total streamfunction into surface and interior components.

        Raises:
            ValueError: If input arrays have invalid shapes or values.
            RuntimeError: If solving the interior streamfunction fails.
        """
        try:
            # Validate inputs
            if self.ssh.ndim != 2 or self.psis.ndim != 3:
                raise ValueError(
                    f"Expected 2D ssh (shape {self.ssh.shape}) and 3D psis (shape {self.psis.shape})"
                )
            if np.any(np.isnan(self.ssh)) or np.any(np.isnan(self.psis)):
                raise ValueError("NaN values detected in ssh or psis")
            ny, nx = self.ssh.shape
            if self.psis.shape[1:] != (ny, nx):
                raise ValueError(
                    f"Shape mismatch: psis shape {self.psis.shape[1:]} vs ssh shape {(ny, nx)}"
                )

            # Compute total streamfunction and anomaly
            self.psit = self.ssh * 9.81 / self.f0
            self.ssha = self.ssh - np.mean(self.ssh)

            # Compute interior streamfunction at the surface
            psii_surface = self.psit - self.psis[-1, :, :]
            if psii_surface.shape != (ny, nx):
                raise ValueError(
                    f"Expected 2D psii_surface, got shape {psii_surface.shape}"
                )

            # Perform QG decomposition
            self.qgdecomposition()

            # Compute coefficients a and b for each grid point
            m1, m2 = 0, 1  # First two eigenmodes
            delta = (
                self.eigenfunction[-1, m1] * self.eigenfunction[0, m2]
                - self.eigenfunction[0, m1] * self.eigenfunction[-1, m2]
            )
            if np.abs(delta) < 1e-10:
                raise ValueError("Singular delta in coefficient computation")

            # Compute coefficients a and b as 2D arrays
            a = (
                psii_surface * self.eigenfunction[0, m2]
                + self.psis[0, :, :] * self.eigenfunction[-1, m2]
            ) / delta
            b = -(
                psii_surface * self.eigenfunction[0, m1]
                + self.psis[0, :, :] * self.eigenfunction[-1, m1]
            ) / delta

            # Vectorized computation of interior streamfunction
            self.psii = np.zeros((self.nz + 3, ny, nx))
            eigen_m1 = self.eigenfunction[:, m1].reshape(-1, 1, 1)
            eigen_m2 = self.eigenfunction[:, m2].reshape(-1, 1, 1)
            self.psii = a * eigen_m1 + b * eigen_m2

            # Compute derived fields
            self.rhoi = self.psi2rho(self.psii)
            self.ui, self.vi = psi2uv(self.lon, self.lat, self.psii)
            self.ut = self.us + self.ui
            self.vt = self.vs + self.vi
            self.vorticity = vorticity(self.lon, self.lat, self.ut, self.vt)

            # Generate RST documentation
            self.rst = [
                f".. _{self.rstfn.replace('.rst', '')}:\n\n",
                "Parameters\n",
                "==================================\n\n",
                f"bottom boundary condition: {self.bottomboundary}\n\n",
                f"size of zc: {self.zc.shape}\n\n",
                f"size of zf: {self.zf.shape}\n\n",
                f"size of n2: {self.n2.shape}\n\n",
                f"size of psis: {self.psis.shape}\n\n",
                f"size of rhos: {self.rhos.shape}\n\n",
                f"size of us,vs: {self.us.shape}\n\n",
                f"size of psii: {self.psii.shape}\n\n",
                f"size of rhoi: {self.rhoi.shape}\n\n",
                f"size of ui,vi: {self.ui.shape}\n\n",
                "================  ================  ================  ================\n",
                " i                  zc                     zf             :math:`N_2` \n",
                "================  ================  ================  ================\n",
            ]
            for i in range(min(self.zc.size, self.zf.size, self.n2.size)):
                self.rst.append(
                    f"{i:16d}{self.zc[i]:16.2f}{self.zf[i]:16.2f}{self.n2[i]:16.2e}\n"
                )
            self.rst.append(
                f"{i+1:16d}{self.zc[i+1]:16.2f}{0:16.2f}{0:16.2e}\n"
                "================  ================  ================  ================\n\n"
            )
            logger.info("Interior streamfunction solution completed")
        except Exception as e:
            logger.error(f"Interior streamfunction solver failed: {e}")
            raise RuntimeError(f"Interior streamfunction error: {e}")

    def plotgrid(self) -> None:
        """
        Plot the vertical grid structure (zc and zf points).

        Creates a schematic of the vertical grid with cell centers (zc) and faces (zf).

        Raises:
            RuntimeError: If plotting fails.
        """
        try:
            fig, ax = plt.subplots()
            ax.plot(np.zeros(4), self.zf[:4], 'ko', markerfacecolor='w', markersize=10)
            for i in range(4):
                ax.text(0.4, self.zf[i] * 1.03, f"$zf_{i}$", fontsize=20)
                ax.text(-0.9, self.zc[i] * 1.02, f"$zc_{i}$", fontsize=20)
            ax.plot(np.zeros(4), self.zc[:4], 'k*')
            ax.plot(
                np.zeros(4), self.zf[-4:] - self.zf[-8], 'ko', markerfacecolor='w', markersize=10
            )
            ax.plot(np.zeros(4), self.zc[-4:] - self.zf[-8], 'k*')
            ax.axis([-5, 5, self.zf[-1] - self.zf[-8] - 50, 50])
            ax.set_yticks([])
            ax.set_xticks([])
            plt.savefig("grid_structure.png")
            plt.close(fig)
            logger.info("Grid plot saved to grid_structure.png")
        except Exception as e:
            logger.error(f"Failed to plot grid: {e}")
            raise RuntimeError(f"Plotting error: {e}")


def fit2Dsurf(x, y, p):
    """
    Given y0=f(t0), find the best fit
    p = a + bx + cy + dx**2 + ey**2 + fxy
    and return a, b, c, d, e, f

    :param x: 1D array of x values
    :param y: 1D array of y values
    :param p: 1D array of z values
    :return: a tuple of coefficients (a, b, c, d, e, f)
    """
    from scipy.optimize import leastsq

    def fit2poly(c, x, y, p):
        a, b, c, d, e, f = c
        err = p - (a + b * x + c * y + d * x**2 + e * y**2 + f * x * y)
        return err

    x = x.flatten()
    y = y.flatten()
    p = p.flatten()
    c = [p.mean(), 1e-3, 1e-3, 1e-6, 1e-6, 1e-6]
    coef = leastsq(fit2poly, c, args=(x, y, p))
    return coef

def fit2exp(x, y, method='exp'):
    """
    Given y0=f(t0), find the best fit p = a + b*exp(c*z) or p = a + b*(1. - tanh(c*z)**2) depending on method,
    and return a,b,c.

    Parameters:
    -----------
    x : array-like
        x-coordinates of the data points
    y : array-like
        y-coordinates of the data points
    method : str, optional
        Type of function to fit. Can be either 'exp' for exponential or 'tanh' for hyperbolic tangent.
        Default is 'exp'.

    Returns:
    --------
    tuple
        Tuple of three values representing the coefficients a, b, and c respectively.

    """
    from scipy.optimize import leastsq
    import numpy as np

    def fit2poly(cc, x, y):
        if method == 'exp':
            err = y - (cc[0] * np.exp(cc[1] * x))
        elif method == 'tanh':
            err = y - (cc[0] * (1. - np.tanh(cc[1] * x) ** 2))
        return err

    x = x.flatten()
    y = y.flatten()
    c = [1e-5, 1. / 200.]
    coef = leastsq(fit2poly, c, args=(x, y))
    return coef[0][0], coef[0][1], coef[0][2]

def anomaly(lon, lat, var):
    """
    Compute the anomaly of a variable with respect to its spatial mean.

    Parameters
    ----------
    lon : array-like
        Array of longitudes.
    lat : array-like
        Array of latitudes.
    var : array-like
        Array of variable values.

    Returns
    -------
    anomaly : array-like
        Anomaly of the variable with respect to its spatial mean.
    """
    from numpy import meshgrid, zeros, arange

    if lon.ndim == 1:
        x1, y1 = meshgrid(lon, lat)
    else:
        x1, y1 = lon, lat

    if var.ndim == 2:
        coef = fit2Dsurf(x1, y1, var)[0]
        varm = fit2poly(coef, x1, y1)
        return var - varm

    elif var.ndim == 3:
        tmp = zeros(var.shape)
        for i in arange(var.shape[0]):
            coef = fit2Dsurf(x1, y1, var[i, :, :])[0]
            varm = fit2poly(coef, x1, y1)
            tmp[i, :, :] = var[i, :, :] - varm

        return tmp
def fit2poly(c,x,y):
    """
    Returns the polynomial fit p = a + bx + cy + dx**2 + ey**2 + fxy
    given the coefficients c and input data x and y.

    Parameters:
        c (list): List of polynomial coefficients
        x (ndarray): Array of x-values
        y (ndarray): Array of y-values

    Returns:
        fit (ndarray): Array of polynomial fit values
    """
    a,b,c,d,e,f=c
    fit = (a + b*x + c*y + d*x**2 + e*y**2 + f*x*y)
    return fit


def lat2f(d):
    """
    Calculates the Coriolis parameter from the latitude d.

    Parameters:
        d (ndarray): Array of latitudes (1- or 2-D)

    Returns:
        f (ndarray): Array of Coriolis parameter values
    """
    from math import pi, sin
    return 2.0*0.729e-4*sin(d*pi/180.0)


def lonlat2xy(lon,lat):
    """
    Converts the given longitude and latitude arrays to x and y coordinates
    using the WGS84 reference ellipsoid.

    Parameters:
        lon (ndarray): Array of longitudes
        lat (ndarray): Array of latitudes

    Returns:
        x (ndarray): Array of x-coordinates
        y (ndarray): Array of y-coordinates
    """
    from pylab import meshgrid,cos,pi
    r = 6371.e3
    lon = lon-lon[0]
    if lon.ndim == 1:
        lon,lat = meshgrid(lon,lat)
    x = 2*pi*r*cos(lat*pi/180.)*lon/360.
    y = 2*pi*r*lat/360.
    return x,y
def vorticity(lon,lat,u,v):
    """
    Calculate vorticity from velocity fields
    lon: longitudes, 1- or 2-D
    lat: latitudes, 1- or 2-D
    u: zonal velocities, same dimensions as lon and lat
    v: meridional velocities, same dimensions as lon and lat
    """
    tmp, dudy = gradxy(lon,lat,u)
    dvdx, tmp = gradxy(lon,lat,v)
    del tmp
    return dvdx-dudy

def gradxy(lon,lat,d):
    """
    Calculate the x and y gradients of a field
    lon: longitudes, 1- or 2-D
    lat: latitudes, 1- or 2-D
    d: field, 2- or 3-D

    The calculated velocity fields share the same dimensions with d,
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
    """
    Calculate horizontal velocity from streamfunction.

    Parameters
    ----------
    lon : 1D or 2D array-like
        Longitudes
    lat : 1D or 2D array-like
        Latitudes
    d : 2D or 3D array-like
        Streamfunction

    Returns
    -------
    u : 2D or 3D array-like
        Zonal velocity
    v : 2D or 3D array-like
        Meridional velocity
    """
    u, v = gradxy(lon, lat, d)
    u = -1 * u
    return u, v

def ssh2uv(lon,lat,ssh):
    """
    Calculate geostrophic velocity from SSH data.

    Parameters:
    -----------
    lon : array-like
        Longitudes.
    lat : array-like
        Latitudes.
    ssh : array-like
        Sea surface height.

    Returns:
    --------
    u : array-like
        Zonal velocity.
    v : array-like
        Meridional velocity.
    """
    from pylab import meshgrid

    if len(lon.shape) == 1:
        lon,lat = meshgrid(lon,lat)

    psi = 9.81 / lat2f(lat) * ssh
    u, v = gradxy(lon,lat,psi)
    u = -1 * u

    return u, v

def twopave(x):
    """
    Compute the two-point average of an array.

    Parameters:
    -----------
    x : array_like
        Input array.

    Returns:
    --------
    y : ndarray
        Output array containing the two-point average of the input array.

    Example:
    --------
    >>> x = [1, 2, 3, 4]
    >>> twopave(x)
    array([1.5, 2.5, 3.5])
    """
    return (x[0:-1]+x[1:])/2.

def N2rho(rhos, n2, dz, rho0=1035., g=9.81, direction='top2bottom'):
    """
    Calculate the mean density field from a surface density field and a stratification profile.

    Parameters:
    -----------
    rhos : array_like
        Surface density field.
    n2 : array_like
        Stratification profile.
    dz : array_like
        Array of thicknesses for each layer.
    rho0 : float, optional
        Reference density, default value is 1035 kg/m^3.
    g : float, optional
        Gravitational acceleration, default value is 9.81 m/s^2.
    direction : str, optional
        Direction of integration, either 'top2bottom' (default) or 'bottom2top'.

    Returns:
    --------
    rho : ndarray
        Mean density field.
    """
    import numpy as np

    nz = n2.size
    ny, nx = rhos.shape
    rho = np.zeros((nz+1, ny, nx))
    rho[0, :, :] = rhos
    dz = abs(dz)

    for i in range(1, nz+1):
        if direction == 'top2bottom':
            rho[i, :, :] = rho[i-1, :, :] + dz[i-1] * n2[i-1] * rho0 / g
        else:
            rho[i, :, :] = rho[i-1, :, :] - dz[i-1] * n2[i-1] * rho0 / g

    return rho
