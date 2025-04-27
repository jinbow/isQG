import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from isqg import isqg_data  

def initialize_isqg_data():
    """Initialize and configure an isqg_data instance."""
    try:
        # Create instance
        data = isqg_data()
        
        # Configuration parameters
        N = 50
        data.z = np.linspace(-500, -10, int(N))
        data.n2 = np.exp(data.z / 200) * 1e-5  # Stratification profile
        data.lon = np.arange(0, 7, 0.5)
        data.lat = np.arange(30, 37, 0.5)
        data.rho0 = 1035.0
        data.filterL = 1e8
        data.bottomboundary = 'rho=0'
        
        return data
    except Exception as e:
        raise RuntimeError(f"Failed to initialize isqg_data: {e}")

def generate_surface_density(data):
    """Generate a Gaussian eddy surface density field."""
    try:
        x, y = data.x, data.y
        x_mean, y_mean = np.mean(x), np.mean(y)
        data.ssd = 0.1 * np.exp(-((x - x_mean) / 1e4) ** 2 - ((y - y_mean) / 1e4) ** 2)
    except Exception as e:
        raise RuntimeError(f"Failed to generate surface density: {e}")

def plot_results(data):
    """Create a 2x2 subplot of the SQG results and a separate grid plot."""
    try:
        # Figure 1: 2x2 subplot
        fig, axs = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
        
        # Plot 1: Surface density
        cf1 = axs[0, 0].contourf(data.ssd, levels=10)
        fig.colorbar(cf1, ax=axs[0, 0])
        axs[0, 0].set_title('Surface Density')
        
        # Plot 2: Density at top level
        cf2 = axs[0, 1].contourf(data.rhos[-1, :, :], levels=10)
        fig.colorbar(cf2, ax=axs[0, 1])
        axs[0, 1].set_title('Density at Top Level')
        
        # Plot 3: Difference between density and surface density
        pc = axs[1, 0].pcolormesh(data.rhos[-1, :, :] - data.ssd)
        fig.colorbar(pc, ax=axs[1, 0])
        axs[1, 0].set_title('Density - Surface Density')
        
        # Plot 4: Vertical cross-section
        cf4 = axs[1, 1].contourf(data.lon, data.zf, data.rhos[:, data.ny // 2, :], levels=10)
        fig.colorbar(cf4, ax=axs[1, 1])
        axs[1, 1].set_title('Vertical Cross-Section')
        
        # Figure 2: Grid plot
        plt.figure(figsize=(6, 6))
        data.plotgrid()
        plt.title('Grid Structure')
        
        plt.show()
    except Exception as e:
        raise RuntimeError(f"Failed to plot results: {e}")

def main():
    """Main function to run the SQG simulation."""
    try:
        # Initialize data
        data = initialize_isqg_data()
        
        # Precondition data
        data.precondition()
        
        # Generate surface density
        generate_surface_density(data)
        
        # Solve SQG equation
        data.solve_sqg()
        
        # Plot results
        plot_results(data)
        
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    main()