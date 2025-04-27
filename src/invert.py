import sys
import logging
from pathlib import Path
from typing import Dict, Any
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from isqg import isqg_data

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_matlab_file(file_path: str) -> Dict[str, Any]:
    """
    Load a MATLAB file and return its contents.

    Args:
        file_path (str): Path to the input MATLAB file.

    Returns:
        Dict[str, Any]: Dictionary containing MATLAB file data.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If the file is not a valid MATLAB file.
    """
    try:
        file_path = Path(file_path)
        if not file_path.is_file():
            raise FileNotFoundError(f"Input file not found: {file_path}")
        if file_path.suffix.lower() != '.mat':
            raise ValueError(f"Input file must be a .mat file, got: {file_path}")
        
        logger.info(f"Loading MATLAB file: {file_path}")
        return loadmat(file_path, squeeze_me=True)
    except Exception as e:
        logger.error(f"Failed to load MATLAB file: {e}")
        raise

def configure_isqg_data(data_dict: Dict[str, Any]) -> isqg_data:
    """
    Configure an isqg_data object with values from a MATLAB data dictionary.

    Args:
        data_dict (Dict[str, Any]): Dictionary containing MATLAB data.

    Returns:
        isqg_data: Configured isqg_data object.

    Raises:
        RuntimeError: If setting attributes fails.
    """
    try:
        data = isqg_data()
        for key, value in data_dict.items():
            if not key.startswith('__'):  # Skip MATLAB metadata
                logger.debug(f"Setting attribute: {key}")
                setattr(data, key, value)
        return data
    except Exception as e:
        logger.error(f"Failed to configure isqg_data: {e}")
        raise RuntimeError(f"Configuration error: {e}")

def solve_equations(data: isqg_data) -> None:
    """
    Solve SQG and PSI equations using the isqg_data object.

    Args:
        data (isqg_data): Configured isqg_data object.

    Raises:
        RuntimeError: If solving equations fails.
    """
    try:
        logger.info("Solving SQG equation")
        data.solve_sqg()
        logger.info("Solving PSI equation")
        data.solve_psii()
    except Exception as e:
        logger.error(f"Failed to solve equations: {e}")
        raise RuntimeError(f"Solving error: {e}")

def plot_diagnostics(data: isqg_data, output_dir: str = ".") -> None:
    """
    Create diagnostic plots for key isqg_data variables.

    Args:
        data (isqg_data): isqg_data object containing results.
        output_dir (str): Directory to save the plot. Defaults to current directory.

    Raises:
        RuntimeError: If plotting fails.
    """
    try:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a 2x2 subplot figure
        fig, axs = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
        axs = axs.flatten()  # Flatten for easier indexing

        # List of variables to plot with their titles and optional vertical cross-section flag
        plot_vars = [
            ("ssd", "Surface Density", False),
            ("rhos", "Density at Top Level", True),  # Plot top level and cross-section
            ("vorticity", "Vorticity", False),
            ("psit", "Total Streamfunction", False)
        ]

        for idx, (var, title, vertical) in enumerate(plot_vars):
            if not hasattr(data, var):
                logger.warning(f"Variable {var} not found in isqg_data object")
                continue

            # Get the variable data
            var_data = getattr(data, var)
            
            # Handle 3D vs 2D data
            if var_data.ndim == 3 and vertical:
                # Plot top level for 3D data
                data_to_plot = var_data[-1, :, :]  # Top level
                if idx < len(axs):
                    cf = axs[idx].contourf(data.lon, data.lat, data_to_plot, levels=10)
                    fig.colorbar(cf, ax=axs[idx])
                    axs[idx].set_title(f"{title} (Top Level)")
                    axs[idx].set_xlabel("Longitude")
                    axs[idx].set_ylabel("Latitude")
                
                # Add vertical cross-section in a separate figure
                fig_v, ax_v = plt.subplots(figsize=(6, 4))
                cf_v = ax_v.contourf(data.lon, data.zf, var_data[:, data.ny // 2, :], levels=10)
                fig_v.colorbar(cf_v, ax=ax_v)
                ax_v.set_title(f"{title} (Vertical Cross-Section)")
                ax_v.set_xlabel("Longitude")
                ax_v.set_ylabel("Depth (m)")
                fig_v.savefig(output_dir / f"{var}_vertical.png")
                plt.close(fig_v)
            elif var_data.ndim == 2 and idx < len(axs):
                # Plot 2D data
                cf = axs[idx].contourf(data.lon, data.lat, var_data, levels=10)
                fig.colorbar(cf, ax=axs[idx])
                axs[idx].set_title(title)
                axs[idx].set_xlabel("Longitude")
                axs[idx].set_ylabel("Latitude")
            else:
                logger.warning(f"Skipping {var}: unexpected data shape {var_data.shape}")

        # Save the main figure
        fig.savefig(output_dir / "diagnostics.png")
        logger.info(f"Diagnostic plots saved to {output_dir / 'diagnostics.png'}")
        plt.close(fig)
        
    except Exception as e:
        logger.error(f"Failed to create diagnostic plots: {e}")
        raise RuntimeError(f"Plotting error: {e}")

def save_results(data: isqg_data, file_path: str) -> None:
    """
    Save specified isqg_data attributes to a MATLAB file.

    Args:
        data (isqg_data): isqg_data object containing results.
        file_path (str): Path to the output MATLAB file.

    Raises:
        RuntimeError: If saving the file fails.
    """
    try:
        file_path = Path(file_path)
        if file_path.suffix.lower() != '.mat':
            raise ValueError(f"Output file must be a .mat file, got: {file_path}")

        output_vars = [
            'us', 'vs', 'ui', 'vi', 'ut', 'vt', 'psis', 'psii', 'psit',
            'ssh', 'ssd', 'rhos', 'rhoi', 'rho0', 'f0', 'z', 'zf', 'zc',
            'lon', 'lat', 'vorticity'
        ]
        
        output_dict = {}
        for var in output_vars:
            if hasattr(data, var):
                output_dict[var] = getattr(data, var)
            else:
                logger.warning(f"Attribute {var} not found in isqg_data object")

        logger.info(f"Saving results to: {file_path}")
        savemat(file_path, output_dict)
    except Exception as e:
        logger.error(f"Failed to save results: {e}")
        raise RuntimeError(f"Save error: {e}")

def invertit(filein: str, fileout: str) -> None:
    """
    Read a MATLAB file with ISQG data, solve SQG and PSI equations, and save results.

    Args:
        filein (str): Path to the input MATLAB file.
        fileout (str): Path to the output MATLAB file.

    Raises:
        RuntimeError: If any step in the process fails.
    """
    try:
        # Load input file
        data_dict = load_matlab_file(filein)
        
        # Configure isqg_data
        data = configure_isqg_data(data_dict)
        
        # Solve equations
        solve_equations(data)
        
        # Plot diagnostics
        output_dir = Path(fileout).parent
        plot_diagnostics(data, output_dir)
        
        # Save results
        save_results(data, fileout)
        
        logger.info("Processing completed successfully")
    except Exception as e:
        logger.error(f"Processing failed: {e}")
        raise RuntimeError(f"Inversion error: {e}")

def main() -> None:
    """
    Command-line interface for the invertit function.

    Expects two command-line arguments: input and output MATLAB file paths.
    """
    try:
        if len(sys.argv) != 3:
            logger.error("Invalid number of arguments")
            print("Usage: python invert.py filein.mat fileout.mat")
            sys.exit(1)
        
        filein, fileout = sys.argv[1:3]
        invertit(filein, fileout)
    except Exception as e:
        logger.error(f"Main execution failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()