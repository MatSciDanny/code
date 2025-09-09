"""
Example usage of the cube visualization script for CRYSTAL calculations.
This shows how to visualize ECH3 (electron density) and POT3 (potential) files.

Run this after you have saved cube_visualizer.py in the same directory.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Import from the cube_visualizer script
try:
    from cube_visualizer import CubeFile, plot_isosurface, plot_slice
except ImportError:
    print("Error: Could not import from cube_visualizer.py")
    print("Make sure cube_visualizer.py is in the same directory as this script")
    sys.exit(1)

# Example 1: Basic electron density visualization
def example_electron_density():
    """Visualize electron charge density from ECH3 calculation."""
    print("Example 1: Electron Density Visualization")
    
    # Read the ECH3 cube file
    cube = CubeFile("material_ECH3.cube")
    cube.print_info()
    
    # Plot with default density isosurfaces
    fig, ax = plot_isosurface(
        cube, 
        iso_values=[0.002, 0.01, 0.05, 0.1],  # e/Bohr³
        cmap='Blues',
        alpha=0.6,
        show_atoms=True
    )
    
    plt.savefig("electron_density_3d.png", dpi=300, bbox_inches='tight')
    plt.show()


# Example 2: Electrostatic potential visualization
def example_potential():
    """Visualize electrostatic potential from POT3 calculation."""
    print("\nExample 2: Electrostatic Potential Visualization")
    
    # Read the POT3 cube file
    cube = CubeFile("material_POT3.cube")
    cube.print_info()
    
    # Plot with positive and negative isosurfaces
    fig, ax = plot_isosurface(
        cube,
        iso_values=[-0.05, -0.01, 0.01, 0.05],  # Hartree
        cmap='RdBu',  # Red for negative, Blue for positive
        alpha=0.7,
        show_atoms=True
    )
    
    plt.savefig("potential_3d.png", dpi=300, bbox_inches='tight')
    plt.show()


# Example 3: 2D slice visualization
def example_2d_slices():
    """Create 2D slices through the volumetric data."""
    print("\nExample 3: 2D Slice Visualization")
    
    cube = CubeFile("material_ECH3.cube")
    
    # Create slices along different directions
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for i, axis in enumerate(['x', 'y', 'z']):
        plt.sca(axes[i])
        
        # Get middle slice position
        n_points = cube.nvoxels[i]
        middle = n_points // 2
        
        # Plot slice
        slice_data = cube.get_slice(axis=axis, position=middle)
        
        # Determine proper extent for the slice
        if axis == 'z':
            extent = [0, cube.nvoxels[0] * cube.voxel_vectors[0, 0] * 0.529177,
                     0, cube.nvoxels[1] * cube.voxel_vectors[1, 1] * 0.529177]
            xlabel, ylabel = 'X (Å)', 'Y (Å)'
        elif axis == 'y':
            extent = [0, cube.nvoxels[0] * cube.voxel_vectors[0, 0] * 0.529177,
                     0, cube.nvoxels[2] * cube.voxel_vectors[2, 2] * 0.529177]
            xlabel, ylabel = 'X (Å)', 'Z (Å)'
        else:  # x
            extent = [0, cube.nvoxels[1] * cube.voxel_vectors[1, 1] * 0.529177,
                     0, cube.nvoxels[2] * cube.voxel_vectors[2, 2] * 0.529177]
            xlabel, ylabel = 'Y (Å)', 'Z (Å)'
        
        im = axes[i].imshow(slice_data.T, origin='lower', extent=extent,
                           cmap='viridis', aspect='auto')
        axes[i].set_xlabel(xlabel)
        axes[i].set_ylabel(ylabel)
        axes[i].set_title(f'Slice along {axis.upper()}')
        
        # Add colorbar
        plt.colorbar(im, ax=axes[i], label='Density (e/Bohr³)')
    
    plt.tight_layout()
    plt.savefig("density_slices.png", dpi=300, bbox_inches='tight')
    plt.show()


# Example 4: Advanced customization
def example_custom_visualization():
    """Advanced visualization with custom settings."""
    print("\nExample 4: Custom Visualization")
    
    # Read cube file
    cube = CubeFile("material_ECH3.cube")
    
    # Create figure with specific size
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get grid data
    X, Y, Z = cube.get_grid_points()
    X *= 0.529177  # Convert to Angstroms
    Y *= 0.529177
    Z *= 0.529177
    
    # Define custom isosurface values based on data statistics
    data_max = np.max(cube.data)
    data_mean = np.mean(cube.data[cube.data > 0])  # Mean of positive values
    
    iso_values = [
        data_mean * 0.1,   # Low density
        data_mean * 0.5,   # Medium density
        data_mean * 2.0,   # High density
        data_max * 0.5     # Very high density
    ]
    
    # Custom colors for each isosurface
    colors = ['lightblue', 'cornflowerblue', 'mediumblue', 'darkblue']
    alphas = [0.3, 0.5, 0.7, 0.9]
    
    # Plot each isosurface
    for iso_val, color, alpha in zip(iso_values, colors, alphas):
        contours = ax.contour3D(X, Y, Z, cube.data, levels=[iso_val],
                               colors=[color], alpha=alpha)
    
    # Plot atoms with custom styling
    element_props = {
        1: {'color': 'white', 'size': 100, 'label': 'H'},    # Hydrogen
        6: {'color': 'black', 'size': 170, 'label': 'C'},    # Carbon
        7: {'color': 'blue', 'size': 155, 'label': 'N'},     # Nitrogen
        8: {'color': 'red', 'size': 152, 'label': 'O'},      # Oxygen
    }
    
    for atomic_num, pos in zip(cube.atomic_numbers, cube.atomic_positions):
        pos_ang = pos * 0.529177  # Convert to Angstroms
        props = element_props.get(atomic_num, {'color': 'gray', 'size': 150, 'label': 'X'})
        
        ax.scatter(pos_ang[0], pos_ang[1], pos_ang[2],
                  c=props['color'], s=props['size'],
                  edgecolors='black', linewidth=2)
    
    # Customize view
    ax.view_init(elev=20, azim=45)
    ax.set_xlabel('X (Å)', fontsize=14)
    ax.set_ylabel('Y (Å)', fontsize=14)
    ax.set_zlabel('Z (Å)', fontsize=14)
    ax.set_title('Custom Electron Density Visualization', fontsize=16)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("custom_visualization.png", dpi=300, bbox_inches='tight')
    plt.show()


# Example 5: Difference density
def example_difference_density():
    """Calculate and visualize difference between two cube files."""
    print("\nExample 5: Difference Density")
    
    # Read two cube files (e.g., different spin states or before/after optimization)
    cube1 = CubeFile("material_ECH3_state1.cube")
    cube2 = CubeFile("material_ECH3_state2.cube")
    
    # Calculate difference
    diff_data = cube2.data - cube1.data
    
    # Create a new cube object for the difference
    cube_diff = cube1  # Copy structure from first cube
    cube_diff.data = diff_data
    cube_diff.data_type = 'difference'
    
    # Plot with positive and negative isosurfaces
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    X, Y, Z = cube_diff.get_grid_points()
    X *= 0.529177
    Y *= 0.529177
    Z *= 0.529177
    
    # Positive difference (density increase) in blue
    pos_levels = [0.001, 0.005, 0.01]
    for level in pos_levels:
        ax.contour3D(X, Y, Z, diff_data, levels=[level],
                    colors=['blue'], alpha=0.6)
    
    # Negative difference (density decrease) in red
    neg_levels = [-0.001, -0.005, -0.01]
    for level in neg_levels:
        ax.contour3D(X, Y, Z, diff_data, levels=[level],
                    colors=['red'], alpha=0.6)
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('Electron Density Difference')
    
    plt.savefig("difference_density.png", dpi=300, bbox_inches='tight')
    plt.show()


# Command line usage examples
def print_usage_examples():
    """Print command line usage examples."""
    print("\n" + "="*60)
    print("COMMAND LINE USAGE EXAMPLES")
    print("="*60 + "\n")
    
    print("1. Basic visualization:")
    print("   python cube_visualizer.py material_ECH3.cube")
    print()
    
    print("2. Custom isosurface values:")
    print("   python cube_visualizer.py material_ECH3.cube --iso 0.001,0.01,0.1")
    print()
    
    print("3. Potential with custom colormap:")
    print("   python cube_visualizer.py material_POT3.cube --cmap RdBu --iso -0.05,0.05")
    print()
    
    print("4. 2D slice visualization:")
    print("   python cube_visualizer.py material_ECH3.cube --slice z 50")
    print()
    
    print("5. Save figure to file:")
    print("   python cube_visualizer.py material_ECH3.cube --save density.png")
    print()
    
    print("6. Hide atoms:")
    print("   python cube_visualizer.py material_ECH3.cube --no-atoms")
    print()
    
    print("7. Transparent visualization:")
    print("   python cube_visualizer.py material_POT3.cube --alpha 0.3")


if __name__ == "__main__":
    # Run examples (comment out those you don't need)
    
    # Uncomment to run specific examples:
    # example_electron_density()
    # example_potential()
    # example_2d_slices()
    # example_custom_visualization()
    # example_difference_density()
    
    # Print usage instructions
    print_usage_examples()
    
    print("\n" + "="*60)
    print("QUICK START")
    print("="*60 + "\n")
    print("To use these scripts:")
    print("1. Save both 'cube_visualizer.py' and this examples file")
    print("2. Place your .cube files in the same directory")
    print("3. Run: python cube_visualizer.py your_file.cube")
    print("4. Or import and use the functions in your own scripts")
    print("\nRequired packages: numpy, matplotlib, scipy")
    print("Optional: scikit-image (for better marching cubes algorithm)")
