import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.io import savemat
from numpy.linalg import cond

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Placeholder for twopave (from ISQGData)
def twopave(arr: np.ndarray) -> np.ndarray:
    """Compute pairwise averages of an array (placeholder)."""
    return (arr[:-1] + arr[1:]) / 2

def construct_matrices(
    z: np.ndarray,
    n2: np.ndarray,
    f0: float,
    plot: bool = False,
    save_file: str = None
) -> dict:
    """
    Construct vertical difference matrices for SQG and eigenvalue decomposition.

    Args:
        z (np.ndarray): Vertical coordinates (m), shape (nz,).
        n2 (np.ndarray): Stratification profile (N^2), shape (nz,), in s^-2.
        f0 (float): Coriolis parameter (s^-1).
        plot (bool): If True, plot grid and matrices. Defaults to False.
        save_file (str): If provided, save results to a MATLAB file. Defaults to None.

    Returns:
        dict: Dictionary containing constructed matrices and arrays:
            - zc: Cell center coordinates (m).
            - zf: Cell face coordinates (m).
            - dzc: Vertical spacing at cell centers (m).
            - dzf: Vertical spacing at cell faces (m).
            - F_eig: Stratification matrix for eigenvalue decomposition.
            - D1_eig: First-order difference matrix for eigenvalue decomposition.
            - D2_eig: Second-order difference matrix for eigenvalue decomposition.
            - eigD2: Eigenvalue operator matrix.
            - F_sqg: Stratification matrix for SQG.
            - D1_sqg: First-order difference matrix for SQG.
            - D2_sqg: Second-order difference matrix for SQG.
            - sqgD2: SQG operator matrix.

    Raises:
        ValueError: If inputs are invalid (e.g., negative n2, non-increasing z).
    """
    try:
        # Validate inputs
        z = np.asarray(z)
        n2 = np.asarray(n2)
        if z.shape != n2.shape or z.ndim != 1:
            raise ValueError(f"Expected 1D z and n2 with same shape, got z={z.shape}, n2={n2.shape}")
        if np.any(n2 <= 0):
            raise ValueError("Non-positive values detected in N^2")
        if np.any(np.isnan(z)) or np.any(np.isnan(n2)):
            raise ValueError("NaN values detected in z or n2")
        if np.any(np.diff(z) <= 0):
            raise ValueError("z must be strictly increasing")
        if f0 == 0:
            raise ValueError("Coriolis parameter (f0) cannot be zero")
        
        nz = len(z)
        logger.info(f"Constructing matrices for nz={nz}")

        # Construct vertical grid
        zf = np.r_[z[0] - (z[1] - z[0]), z, z[-1] + (z[-1] - z[-2])]  # Shape (nz+2,)
        zc = np.r_[
            zf[0] - 0.5 * (zf[1] - zf[0]),
            twopave(zf),
            zf[-1] + 0.5 * (zf[-1] - zf[-2])
        ]  # Shape (nz+3,)
        dzc = np.diff(zc)  # Shape (nz+2,)
        dzf = np.diff(zf)  # Shape (nz+1,)
        
        logger.debug(f"Shapes: zf={zf.shape}, zc={zc.shape}, dzc={dzc.shape}, dzf={dzf.shape}")

        # Eigenvalue decomposition matrices
        F_eig_values = np.r_[0, f0**2 / n2, 0]  # Shape (nz+2,)
        F_eig = sparse.diags(F_eig_values, format='csr')  # Shape (nz+2, nz+2)
        D1_eig = sparse.diags(1. / dzc, format='csr') @ sparse.csr_matrix(
            np.diff(-np.identity(nz + 2), axis=1)
        )  # Shape (nz+2, nz+1)
        D2_eig = sparse.diags(1. / dzf, format='csr') @ sparse.csr_matrix(
            np.diff(np.identity(nz + 2), axis=0)
        )  # Shape (nz+1, nz+2)
        eigD2 = D2_eig @ F_eig @ D1_eig  # Shape (nz+1, nz+1)
        
        logger.debug(f"Eigen matrices: F_eig={F_eig.shape}, D1_eig={D1_eig.shape}, D2_eig={D2_eig.shape}, eigD2={eigD2.shape}")

        # SQG matrices
        F_sqg_values = np.r_[f0**2 / n2[0], f0**2 / n2, f0**2 / n2[-1]]  # Shape (nz+3,)
        F_sqg = sparse.diags(F_sqg_values, format='csr')  # Shape (nz+3, nz+3)
        D1_sqg = sparse.diags(1. / dzc, format='csr') @ sparse.csr_matrix(
            np.diff(-np.identity(nz + 3), axis=0)
        )  # Shape (nz+2, nz+3)
        D2_sqg = sparse.diags(np.r_[1., 1. / dzf, 1.], format='csr') @ sparse.csr_matrix(
            np.diff(np.identity(nz + 3), axis=1)
        )  # Shape (nz+3, nz+2)
        sqgD2 = D2_sqg @ F_sqg @ D1_sqg  # Shape (nz+3, nz+3)
        
        # Apply boundary conditions
        sqgD2_lil = sqgD2.tolil()
        sqgD2_lil[-1, -2:] = [-1. / dzc[-1], 1. / dzc[-1]]
        sqgD2_lil[0, :2] = [-1. / dzc[0], 1. / dzc[0]]
        sqgD2 = sqgD2_lil.tocsr()
        
        logger.debug(f"SQG matrices: F_sqg={F_sqg.shape}, D1_sqg={D1_sqg.shape}, D2_sqg={D2_sqg.shape}, sqgD2={sqgD2.shape}")

        # Check numerical stability
        eigD2_dense = eigD2.toarray()
        sqgD2_dense = sqgD2.toarray()
        if cond(eigD2_dense) > 1e12:
            logger.warning(f"Ill-conditioned eigD2 matrix, condition number: {cond(eigD2_dense)}")
        if cond(sqgD2_dense) > 1e12:
            logger.warning(f"Ill-conditioned sqgD2 matrix, condition number: {cond(sqgD2_dense)}")

        # Store results
        results = {
            'zc': zc,
            'zf': zf,
            'dzc': dzc,
            'dzf': dzf,
            'F_eig': F_eig,
            'D1_eig': D1_eig,
            'D2_eig': D2_eig,
            'eigD2': eigD2,
            'F_sqg': F_sqg,
            'D1_sqg': D1_sqg,
            'D2_sqg': D2_sqg,
            'sqgD2': sqgD2
        }

        # Plotting (optional)
        if plot:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            
            # Plot grid
            axes[0].plot(np.zeros_like(zf), zf, 'ko', label='zf')
            axes[0].plot(np.zeros_like(zc), zc, 'r*', label='zc')
            axes[0].set_xlabel('Index')
            axes[0].set_ylabel('Depth (m)')
            axes[0].set_title('Vertical Grid')
            axes[0].legend()
            axes[0].grid(True)
            
            # Plot eigD2
            im1 = axes[1].imshow(eigD2_dense, cmap='viridis')
            axes[1].set_title('eigD2 Matrix')
            plt.colorbar(im1, ax=axes[1])
            
            # Plot sqgD2
            im2 = axes[2].imshow(sqgD2_dense, cmap='viridis')
            axes[2].set_title('sqgD2 Matrix')
            plt.colorbar(im2, ax=axes[2])
            
            plt.tight_layout()
            plt.savefig('matrix_visualization.png')
            plt.close(fig)
            logger.info("Matrix visualization saved to matrix_visualization.png")

        # Save to MATLAB file (optional)
        if save_file:
            savemat(save_file, {k: v.toarray() if sparse.issparse(v) else v for k, v in results.items()})
            logger.info(f"Results saved to {save_file}")

        return results

    except Exception as e:
        logger.error(f"Matrix construction failed: {e}")
        raise RuntimeError(f"Matrix construction error: {e}")

def main():
    """Main function to test matrix construction."""
    # Define test parameters
    z = np.arange(-100, -10, 10)  # Depth levels (m)
    n2 = np.ones_like(z) * 1.0   # Constant N^2 (s^-2)
    f0 = 1.0                     # Coriolis parameter (s^-1)

    # Construct matrices
    results = construct_matrices(z, n2, f0, plot=True, save_file='test_matrices.mat')

    # Log results
    logger.info("Matrix construction completed successfully")
    for key, value in results.items():
        shape = value.shape if hasattr(value, 'shape') else 'sparse'
        logger.debug(f"{key}: shape={shape}")

if __name__ == "__main__":
    main()