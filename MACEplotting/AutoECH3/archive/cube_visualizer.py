"""
Script to read and visualize cube files from CRYSTAL calculations.
Supports ECH3 (electron density) and POT3 (electrostatic potential) files.

Usage:
python cube_visualizer.py filename.cube [options]

Options:
--iso VALUES        Isosurface values (comma-separated)
--cmap COLORMAP     Matplotlib colormap name (default: 'RdBu')
--alpha ALPHA       Transparency (0-1, default: 0.7)
--slice AXIS POS    Show slice along AXIS at position POS
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import ndimage
import sys
import os
import argparse

class CubeFile:
    """Class to handle cube file reading and data storage."""
    
    def __init__(self, filename):
        self.filename = filename
        self.comment1 = ""
        self.comment2 = ""
        self.natoms = 0
        self.origin = np.zeros(3)
        self.nvoxels = np.zeros(3, dtype=int)
        self.voxel_vectors = np.zeros((3, 3))
        self.atoms = []
        self.atomic_numbers = []
        self.atomic_positions = []
        self.data = None
        self.data_type = None  # 'density' or 'potential'
        
        self.read_cube()
        self._determine_data_type()
    
    def read_cube(self):
        """Read cube file and extract all information."""
        with open(self.filename, 'r') as f:
            # Read header
            self.comment1 = f.readline().strip()
            self.comment2 = f.readline().strip()
            
            # Read atoms and origin
            line = f.readline().split()
            self.natoms = int(line[0])
            self.origin = np.array([float(line[1]), float(line[2]), float(line[3])])
            
            # Read voxel information
            for i in range(3):
                line = f.readline().split()
                self.nvoxels[i] = int(line[0])
                self.voxel_vectors[i] = np.array([float(line[1]), float(line[2]), float(line[3])])
            
            # Read atomic positions
            for i in range(abs(self.natoms)):
                line = f.readline().split()
                atomic_number = int(line[0])
                charge = float(line[1])
                position = np.array([float(line[2]), float(line[3]), float(line[4])])
                self.atomic_numbers.append(atomic_number)
                self.atomic_positions.append(position)
            
            # Read volumetric data
            data_list = []
            for line in f:
                data_list.extend([float(x) for x in line.split()])
            
            # Reshape data to 3D array
            self.data = np.array(data_list).reshape(self.nvoxels)
    
    def _determine_data_type(self):
        """Determine if data is electron density or potential based on filename and values."""
        filename_lower = self.filename.lower()
        if 'ech3' in filename_lower or 'density' in filename_lower:
            self.data_type = 'density'
        elif 'pot3' in filename_lower or 'potential' in filename_lower:
            self.data_type = 'potential'
        else:
            # Try to guess from data characteristics
            if np.all(self.data >= 0):
                self.data_type = 'density'
            else:
                self.data_type = 'potential'
    
    def get_grid_points(self):
        """Generate 3D grid points for the cube data."""
        x = np.arange(self.nvoxels[0]) * self.voxel_vectors[0, 0] + self.origin[0]
        y = np.arange(self.nvoxels[1]) * self.voxel_vectors[1, 1] + self.origin[1]
        z = np.arange(self.nvoxels[2]) * self.voxel_vectors[2, 2] + self.origin[2]
        return np.meshgrid(x, y, z, indexing='ij')
    
    def get_slice(self, axis='z', position=None):
        """Extract a 2D slice from the 3D data."""
        if position is None:
            # Default to middle slice
            if axis == 'x':
                position = self.nvoxels[0] // 2
            elif axis == 'y':
                position = self.nvoxels[1] // 2
            else:  # z
                position = self.nvoxels[2] // 2
        
        if axis == 'x':
            slice_data = self.data[position, :, :]
        elif axis == 'y':
            slice_data = self.data[:, position, :]
        else:  # z
            slice_data = self.data[:, :, position]
        
        return slice_data
    
    def print_info(self):
        """Print summary information about the cube file."""
        print(f"\nCube File: {self.filename}")
        print(f"Comments: {self.comment1} | {self.comment2}")
        print(f"Data type: {self.data_type}")
        print(f"Number of atoms: {abs(self.natoms)}")
        print(f"Origin: {self.origin} Bohr")
        print(f"Grid points: {self.nvoxels[0]} x {self.nvoxels[1]} x {self.nvoxels[2]}")
        print(f"Voxel dimensions: {self.voxel_vectors[0, 0]:.4f} x {self.voxel_vectors[1, 1]:.4f} x {self.voxel_vectors[2, 2]:.4f} Bohr")
        print(f"Data range: [{np.min(self.data):.6f}, {np.max(self.data):.6f}]")
        print(f"Data mean: {np.mean(self.data):.6f}")
        print(f"Data std: {np.std(self.data):.6f}")


def plot_isosurface(cube, iso_values=None, cmap='RdBu', alpha=0.7, show_atoms=True):
    """Plot 3D isosurfaces of cube data."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Get grid points
    X, Y, Z = cube.get_grid_points()
    
    # Convert from Bohr to Angstrom for display
    bohr_to_ang = 0.529177
    X *= bohr_to_ang
    Y *= bohr_to_ang
    Z *= bohr_to_ang
    
    # Determine isosurface values if not provided
    if iso_values is None:
        if cube.data_type == 'density':
            # For density, use log-spaced values
            max_val = np.max(cube.data)
            iso_values = [0.001, 0.01, 0.05]
        else:
            # For potential, use symmetric positive/negative values
            max_abs = np.max(np.abs(cube.data))
            iso_values = [-0.1, -0.01, 0.01, 0.1]
    
    # Create colormap
    colors = plt.cm.get_cmap(cmap)
    norm = plt.Normalize(vmin=min(iso_values), vmax=max(iso_values))
    
    # Plot isosurfaces
    for iso_val in iso_values:
        if cube.data_type == 'density' and iso_val < 0:
            continue  # Skip negative values for density
        
        try:
            # Use marching cubes algorithm (simplified version)
            from skimage import measure
            verts, faces, normals, values = measure.marching_cubes(
                cube.data, level=iso_val, spacing=(
                    cube.voxel_vectors[0, 0] * bohr_to_ang,
                    cube.voxel_vectors[1, 1] * bohr_to_ang,
                    cube.voxel_vectors[2, 2] * bohr_to_ang
                )
            )
            
            # Adjust vertices to account for origin
            verts[:, 0] += cube.origin[0] * bohr_to_ang
            verts[:, 1] += cube.origin[1] * bohr_to_ang
            verts[:, 2] += cube.origin[2] * bohr_to_ang
            
            # Plot surface
            color = colors(norm(iso_val))
            ax.plot_trisurf(
                verts[:, 0], verts[:, 1], faces, verts[:, 2],
                color=color, alpha=alpha, shade=True,
                label=f'Iso = {iso_val:.3f}'
            )
        except:
            # Fallback to contour3d if skimage is not available
            contours = ax.contour3D(X, Y, Z, cube.data, levels=[iso_val],
                                   colors=[colors(norm(iso_val))], alpha=alpha)
    
    # Plot atoms if requested
    if show_atoms:
        atom_colors = {
            1: 'white',    # H
            6: 'black',    # C
            7: 'blue',     # N
            8: 'red',      # O
            9: 'green',    # F
            15: 'orange',  # P
            16: 'yellow',  # S
        }
        
        for i, (atomic_num, pos) in enumerate(zip(cube.atomic_numbers, cube.atomic_positions)):
            pos_ang = pos * bohr_to_ang
            color = atom_colors.get(atomic_num, 'gray')
            size = 100 * (atomic_num / 10 + 1)  # Scale size by atomic number
            ax.scatter(pos_ang[0], pos_ang[1], pos_ang[2], 
                      c=color, s=size, edgecolors='black', linewidth=1)
    
    # Set labels and title
    ax.set_xlabel('X (Å)', fontsize=12)
    ax.set_ylabel('Y (Å)', fontsize=12)
    ax.set_zlabel('Z (Å)', fontsize=12)
    
    title = f"{cube.data_type.capitalize()} - {os.path.basename(cube.filename)}"
    ax.set_title(title, fontsize=14)
    
    # Set aspect ratio
    ax.set_box_aspect([1,1,1])
    
    # Add legend
    ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    return fig, ax


def plot_slice(cube, axis='z', position=None, cmap='RdBu', show_contours=True):
    """Plot 2D slice of cube data."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Get slice data
    slice_data = cube.get_slice(axis, position)
    
    # Determine extent based on axis
    bohr_to_ang = 0.529177
    if axis == 'z':
        extent = [
            cube.origin[0] * bohr_to_ang,
            (cube.origin[0] + cube.nvoxels[0] * cube.voxel_vectors[0, 0]) * bohr_to_ang,
            cube.origin[1] * bohr_to_ang,
            (cube.origin[1] + cube.nvoxels[1] * cube.voxel_vectors[1, 1]) * bohr_to_ang
        ]
        xlabel, ylabel = 'X (Å)', 'Y (Å)'
    elif axis == 'y':
        extent = [
            cube.origin[0] * bohr_to_ang,
            (cube.origin[0] + cube.nvoxels[0] * cube.voxel_vectors[0, 0]) * bohr_to_ang,
            cube.origin[2] * bohr_to_ang,
            (cube.origin[2] + cube.nvoxels[2] * cube.voxel_vectors[2, 2]) * bohr_to_ang
        ]
        xlabel, ylabel = 'X (Å)', 'Z (Å)'
    else:  # x
        extent = [
            cube.origin[1] * bohr_to_ang,
            (cube.origin[1] + cube.nvoxels[1] * cube.voxel_vectors[1, 1]) * bohr_to_ang,
            cube.origin[2] * bohr_to_ang,
            (cube.origin[2] + cube.nvoxels[2] * cube.voxel_vectors[2, 2]) * bohr_to_ang
        ]
        xlabel, ylabel = 'Y (Å)', 'Z (Å)'
    
    # Plot image
    im = ax.imshow(slice_data.T, origin='lower', extent=extent, cmap=cmap, aspect='auto')
    
    # Add contours if requested
    if show_contours:
        X = np.linspace(extent[0], extent[1], slice_data.shape[0])
        Y = np.linspace(extent[2], extent[3], slice_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        
        if cube.data_type == 'density':
            levels = np.logspace(-3, 0, 10)
            colors = 'black'
        else:
            levels = 15
            colors = 'black'
        
        contours = ax.contour(X, Y, slice_data.T, levels=levels, colors=colors, 
                             alpha=0.5, linewidths=0.5)
        ax.clabel(contours, inline=True, fontsize=8, fmt='%.3f')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label=f'{cube.data_type.capitalize()}')
    
    # Set labels and title
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    
    position_ang = position * cube.voxel_vectors[['x', 'y', 'z'].index(axis)][['x', 'y', 'z'].index(axis)] * bohr_to_ang
    title = f"{cube.data_type.capitalize()} - Slice along {axis.upper()} at {position_ang:.2f} Å"
    ax.set_title(title, fontsize=14)
    
    plt.tight_layout()
    return fig, ax


def main():
    """Main function to handle command line arguments and plotting."""
    parser = argparse.ArgumentParser(description='Visualize cube files from CRYSTAL calculations')
    parser.add_argument('filename', help='Cube file to visualize')
    parser.add_argument('--iso', type=str, help='Comma-separated isosurface values')
    parser.add_argument('--cmap', default='RdBu', help='Colormap name')
    parser.add_argument('--alpha', type=float, default=0.7, help='Transparency (0-1)')
    parser.add_argument('--slice', nargs=2, help='Slice axis and position')
    parser.add_argument('--no-atoms', action='store_true', help='Do not show atoms')
    parser.add_argument('--save', help='Save figure to file')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.filename):
        print(f"Error: File '{args.filename}' not found!")
        sys.exit(1)
    
    # Read cube file
    print(f"Reading cube file: {args.filename}")
    cube = CubeFile(args.filename)
    cube.print_info()
    
    # Parse isosurface values
    iso_values = None
    if args.iso:
        iso_values = [float(x) for x in args.iso.split(',')]
    
    # Create plots
    if args.slice:
        # 2D slice plot
        axis = args.slice[0].lower()
        position = int(args.slice[1])
        fig, ax = plot_slice(cube, axis=axis, position=position, cmap=args.cmap)
    else:
        # 3D isosurface plot
        fig, ax = plot_isosurface(cube, iso_values=iso_values, cmap=args.cmap, 
                                 alpha=args.alpha, show_atoms=not args.no_atoms)
    
    # Save or show figure
    if args.save:
        plt.savefig(args.save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {args.save}")
    else:
        plt.show()


# Example usage functions for interactive use
def quick_plot_density(filename, iso_values=[0.001, 0.01, 0.1]):
    """Quick function to plot electron density."""
    cube = CubeFile(filename)
    cube.print_info()
    fig, ax = plot_isosurface(cube, iso_values=iso_values, cmap='Blues', alpha=0.7)
    plt.show()
    return cube, fig, ax


def quick_plot_potential(filename, iso_values=[-0.1, -0.01, 0.01, 0.1]):
    """Quick function to plot electrostatic potential."""
    cube = CubeFile(filename)
    cube.print_info()
    fig, ax = plot_isosurface(cube, iso_values=iso_values, cmap='RdBu', alpha=0.7)
    plt.show()
    return cube, fig, ax


def compare_cube_files(file1, file2, iso_value=0.01):
    """Compare two cube files side by side."""
    fig = plt.figure(figsize=(16, 8))
    
    # First cube
    ax1 = fig.add_subplot(121, projection='3d')
    cube1 = CubeFile(file1)
    X1, Y1, Z1 = cube1.get_grid_points()
    bohr_to_ang = 0.529177
    X1 *= bohr_to_ang
    Y1 *= bohr_to_ang  
    Z1 *= bohr_to_ang
    ax1.contour3D(X1, Y1, Z1, cube1.data, levels=[iso_value], colors='blue', alpha=0.7)
    ax1.set_title(f"{os.path.basename(file1)}")
    
    # Second cube
    ax2 = fig.add_subplot(122, projection='3d')
    cube2 = CubeFile(file2)
    X2, Y2, Z2 = cube2.get_grid_points()
    X2 *= bohr_to_ang
    Y2 *= bohr_to_ang
    Z2 *= bohr_to_ang
    ax2.contour3D(X2, Y2, Z2, cube2.data, levels=[iso_value], colors='red', alpha=0.7)
    ax2.set_title(f"{os.path.basename(file2)}")
    
    plt.tight_layout()
    plt.show()
    return cube1, cube2, fig


if __name__ == "__main__":
    main()
