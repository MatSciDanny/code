#!/usr/bin/env python3
"""
Complete CRYSTAL Cube File Visualizer
Supports ECH3 (electron density), POT3 (electrostatic potential), and SPIN (spin density) files.
Automatically detects 2D vs 3D structures and suggests appropriate visualization parameters.

Usage:
python crystal_cube_viz.py filename.cube [options]

Options:
--iso VALUES        Override isosurface values (comma-separated)
--cmap COLORMAP     Matplotlib colormap name (default: auto-detected)
--alpha ALPHA       Transparency (0-1, default: 0.7)
--slice AXIS POS    Show 2D slice along AXIS at position POS
--no-atoms         Don't show atoms
--save FILENAME    Save figure to file
--view ANGLES      Set viewing angles (elev,azim) for 3D plot
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
import os
import argparse
from typing import Tuple, List, Optional, Dict

# Try to import skimage for better isosurface rendering
try:
    from skimage import measure
    HAS_SKIMAGE = True
except ImportError:
    HAS_SKIMAGE = False
    print("Note: Install scikit-image for better isosurface rendering")
    print("      pip install scikit-image")
    print()


class CubeFile:
    """Class to handle cube file reading and data storage."""
    
    def __init__(self, filename: str):
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
        self.data_type = None
        self.is_2d = False
        self.z_vacuum = 0.0
        self.cell_parameters = {}
        
        self.read_cube()
        self._determine_data_type()
        self._analyze_structure()
        self._get_recommended_isovalues()
    
    def read_cube(self):
        """Read cube file and extract all information."""
        with open(self.filename, 'r') as f:
            # Read header
            self.comment1 = f.readline().strip()
            self.comment2 = f.readline().strip()
            
            # Extract cell parameters from comment if available
            self._parse_cell_parameters()
            
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
    
    def _parse_cell_parameters(self):
        """Extract cell parameters from comment line."""
        # Look for cell parameters in format: a b c alpha beta gamma
        try:
            parts = self.comment2.split()
            if len(parts) >= 6:
                # Try to find 6 consecutive numbers
                for i in range(len(parts) - 5):
                    try:
                        a = float(parts[i])
                        b = float(parts[i+1])
                        c = float(parts[i+2])
                        alpha = float(parts[i+3])
                        beta = float(parts[i+4])
                        gamma = float(parts[i+5])
                        
                        # Sanity check - cell parameters should be reasonable
                        if (0 < a < 1000 and 0 < b < 1000 and 0 < c < 10000 and
                            0 < alpha <= 180 and 0 < beta <= 180 and 0 < gamma <= 180):
                            self.cell_parameters = {
                                'a': a, 'b': b, 'c': c,
                                'alpha': alpha, 'beta': beta, 'gamma': gamma
                            }
                            break
                    except ValueError:
                        continue
        except:
            pass
    
    def _determine_data_type(self):
        """Determine if data is electron density, spin density, or potential."""
        filename_lower = self.filename.lower()
        if 'spin' in filename_lower:
            self.data_type = 'spin'
        elif 'ech3' in filename_lower or ('density' in filename_lower and 'spin' not in filename_lower):
            self.data_type = 'density'
        elif 'pot' in filename_lower or 'potential' in filename_lower:
            self.data_type = 'potential'
        else:
            # Guess from data characteristics
            if np.all(self.data >= 0):
                self.data_type = 'density'
            else:
                self.data_type = 'spin' if np.abs(np.mean(self.data)) < 1e-6 else 'potential'
    
    def _analyze_structure(self):
        """Analyze if structure is 2D or 3D based on cell and atomic positions."""
        # Get cell dimensions
        cell_x = self.nvoxels[0] * np.linalg.norm(self.voxel_vectors[0])
        cell_y = self.nvoxels[1] * np.linalg.norm(self.voxel_vectors[1])
        cell_z = self.nvoxels[2] * np.linalg.norm(self.voxel_vectors[2])
        
        # Use cell parameters if available
        if self.cell_parameters:
            c_param = self.cell_parameters.get('c', 0)
            # If c parameter is very large (>20 Bohr), likely 2D
            if c_param > 20:
                self.is_2d = True
                self.z_vacuum = c_param
        
        # Check z-extension of atoms
        z_positions = [pos[2] for pos in self.atomic_positions]
        z_min, z_max = min(z_positions), max(z_positions)
        z_extent = z_max - z_min
        
        # If z-extent is small compared to cell_z, likely 2D
        if z_extent < 5.0 and cell_z > 20.0:  # in Bohr
            self.is_2d = True
            self.z_vacuum = cell_z
        
        # For 2D materials, adjust atomic positions if needed
        if self.is_2d:
            z_center = (z_min + z_max) / 2
            # If atoms are split by periodic boundary, bring them together
            if z_extent > cell_z / 2:
                # Atoms are split across boundary - need to wrap
                new_positions = []
                for pos in self.atomic_positions:
                    new_pos = pos.copy()
                    if pos[2] < z_center - cell_z / 4:
                        new_pos[2] += cell_z
                    new_positions.append(new_pos)
                self.atomic_positions = new_positions
    
    def _get_recommended_isovalues(self) -> List[float]:
        """Calculate recommended isovalues based on data statistics."""
        # Exclude zeros for better statistics
        nonzero_data = self.data[self.data != 0]
        
        if len(nonzero_data) == 0:
            # All zeros - use very small values
            self.recommended_iso = {
                'density': [1e-6, 1e-5, 1e-4],
                'potential': [-1e-4, -1e-5, 1e-5, 1e-4],
                'spin': [-1e-5, -1e-6, 1e-6, 1e-5]
            }.get(self.data_type, [1e-5])
            return
        
        # Calculate statistics
        data_min = np.min(self.data)
        data_max = np.max(self.data)
        data_mean = np.mean(nonzero_data)
        data_std = np.std(nonzero_data)
        max_abs = np.max(np.abs(self.data))
        
        if self.data_type == 'density':
            # Log-spaced values for density
            if data_max > 0.1:
                self.recommended_iso = [0.001, 0.01, 0.05, 0.1]
            elif data_max > 0.01:
                self.recommended_iso = [0.0001, 0.001, 0.01, 0.05]
            else:
                self.recommended_iso = [1e-5, 1e-4, 0.001, 0.01]
        
        elif self.data_type == 'spin':
            # Symmetric values for spin density
            if max_abs < 1e-6:
                self.recommended_iso = [-1e-7, -1e-8, 1e-8, 1e-7]
            elif max_abs < 1e-3:
                self.recommended_iso = [-1e-4, -1e-5, 1e-5, 1e-4]
            elif max_abs < 0.01:
                self.recommended_iso = [-0.001, -0.0001, 0.0001, 0.001]
            else:
                self.recommended_iso = [-0.01, -0.001, 0.001, 0.01]
        
        else:  # potential
            # Symmetric values around zero
            levels = sorted([max_abs * 0.1, max_abs * 0.3, max_abs * 0.5])
            self.recommended_iso = []
            for level in levels:
                if level > 1e-6:
                    self.recommended_iso.extend([-level, level])
            self.recommended_iso = sorted(self.recommended_iso)[:4]
    
    def get_grid_points(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate 3D grid points for the cube data."""
        x = np.arange(self.nvoxels[0]) * self.voxel_vectors[0, 0] + self.origin[0]
        y = np.arange(self.nvoxels[1]) * self.voxel_vectors[1, 1] + self.origin[1]
        z = np.arange(self.nvoxels[2]) * self.voxel_vectors[2, 2] + self.origin[2]
        return np.meshgrid(x, y, z, indexing='ij')
    
    def get_slice(self, axis: str = 'z', position: Optional[int] = None) -> np.ndarray:
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
        print(f"Data type: {self.data_type}")
        print(f"Structure type: {'2D' if self.is_2d else '3D'}")
        print(f"Number of atoms: {abs(self.natoms)}")
        print(f"Grid points: {self.nvoxels[0]} × {self.nvoxels[1]} × {self.nvoxels[2]}")
        
        if self.cell_parameters:
            print(f"Cell parameters: a={self.cell_parameters['a']:.3f} b={self.cell_parameters['b']:.3f} "
                  f"c={self.cell_parameters['c']:.3f} Bohr")
            print(f"Angles: α={self.cell_parameters['alpha']:.1f}° β={self.cell_parameters['beta']:.1f}° "
                  f"γ={self.cell_parameters['gamma']:.1f}°")
        
        nonzero = self.data[self.data != 0]
        if len(nonzero) > 0:
            print(f"Data range: [{np.min(self.data):.6f}, {np.max(self.data):.6f}]")
            print(f"Mean (nonzero): {np.mean(nonzero):.6f}")
            print(f"Recommended isovalues: {', '.join([f'{x:.2e}' for x in self.recommended_iso])}")
        else:
            print("Data contains all zeros!")


def plot_isosurface(cube: CubeFile, iso_values: Optional[List[float]] = None, 
                   cmap: Optional[str] = None, alpha: float = 0.7, 
                   show_atoms: bool = True, view_angles: Optional[Tuple[float, float]] = None,
                   fig_ax: Optional[Tuple] = None) -> Tuple:
    """Plot 3D isosurfaces of cube data."""
    
    # Create figure if not provided
    if fig_ax is None:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig, ax = fig_ax
    
    # Get grid points
    X, Y, Z = cube.get_grid_points()
    
    # Convert from Bohr to Angstrom
    bohr_to_ang = 0.529177
    X *= bohr_to_ang
    Y *= bohr_to_ang
    Z *= bohr_to_ang
    
    # Use recommended values if not provided
    if iso_values is None:
        iso_values = cube.recommended_iso
        print(f"Using recommended isovalues: {', '.join([f'{x:.2e}' for x in iso_values])}")
    
    # Choose colormap
    if cmap is None:
        cmap = {
            'density': 'Blues',
            'spin': 'RdBu',
            'potential': 'RdBu'
        }.get(cube.data_type, 'viridis')
    
    # Create colormap
    colors_map = plt.cm.get_cmap(cmap)
    norm = plt.Normalize(vmin=min(iso_values), vmax=max(iso_values))
    
    # Plot isosurfaces
    plotted_surfaces = 0
    for iso_val in iso_values:
        if cube.data_type == 'density' and iso_val < 0:
            continue
        
        # Check if this isosurface exists
        if not np.any((cube.data >= iso_val) if iso_val > 0 else (cube.data <= iso_val)):
            continue
        
        color = colors_map(norm(iso_val))
        
        if HAS_SKIMAGE and plotted_surfaces < 5:  # Limit for performance
            try:
                verts, faces, normals, values = measure.marching_cubes(
                    cube.data, level=iso_val, spacing=(
                        cube.voxel_vectors[0, 0] * bohr_to_ang,
                        cube.voxel_vectors[1, 1] * bohr_to_ang,
                        cube.voxel_vectors[2, 2] * bohr_to_ang
                    )
                )
                
                verts[:, 0] += cube.origin[0] * bohr_to_ang
                verts[:, 1] += cube.origin[1] * bohr_to_ang
                verts[:, 2] += cube.origin[2] * bohr_to_ang
                
                ax.plot_trisurf(
                    verts[:, 0], verts[:, 1], faces, verts[:, 2],
                    color=color, alpha=alpha, shade=True
                )
                plotted_surfaces += 1
            except:
                ax.contour3D(X, Y, Z, cube.data, levels=[iso_val],
                           colors=[color], alpha=alpha)
        else:
            ax.contour3D(X, Y, Z, cube.data, levels=[iso_val],
                       colors=[color], alpha=alpha)
    
    # Plot atoms
    if show_atoms:
        atom_colors = {
            1: ('white', 'H'), 6: ('black', 'C'), 7: ('blue', 'N'), 
            8: ('red', 'O'), 9: ('green', 'F'), 14: ('tan', 'Si'),
            15: ('orange', 'P'), 16: ('yellow', 'S'), 17: ('green', 'Cl'),
            34: ('orange', 'Se'), 50: ('gray', 'Sn'), 52: ('darkgoldenrod', 'Te')
        }
        
        for atomic_num, pos in zip(cube.atomic_numbers, cube.atomic_positions):
            pos_ang = pos * bohr_to_ang
            
            # For 2D materials, ensure atoms are shown together
            if cube.is_2d and pos_ang[2] > cube.z_vacuum * bohr_to_ang / 2:
                pos_ang[2] -= cube.z_vacuum * bohr_to_ang
            
            color, element = atom_colors.get(atomic_num, ('gray', 'X'))
            size = 100 * (1 + atomic_num / 20)
            ax.scatter(pos_ang[0], pos_ang[1], pos_ang[2], 
                      c=color, s=size, edgecolors='black', linewidth=1)
    
    # Set labels and title
    ax.set_xlabel('X (Å)', fontsize=12)
    ax.set_ylabel('Y (Å)', fontsize=12)
    ax.set_zlabel('Z (Å)', fontsize=12)
    
    title_parts = [cube.data_type.capitalize()]
    if cube.is_2d:
        title_parts.append("(2D material)")
    title_parts.append(f"- {os.path.basename(cube.filename)}")
    ax.set_title(' '.join(title_parts), fontsize=14)
    
    # Set aspect ratio - for 2D materials, compress Z
    if cube.is_2d:
        # Calculate appropriate Z limits to show just the material
        z_positions = [pos[2] * bohr_to_ang for pos in cube.atomic_positions]
        z_min, z_max = min(z_positions), max(z_positions)
        z_center = (z_min + z_max) / 2
        z_range = max(z_max - z_min, 5.0)  # At least 5 Å range
        ax.set_zlim(z_center - z_range, z_center + z_range)
    else:
        ax.set_box_aspect([1, 1, 1])
    
    # Set view angle
    if view_angles:
        ax.view_init(elev=view_angles[0], azim=view_angles[1])
    elif cube.is_2d:
        ax.view_init(elev=20, azim=45)  # Good angle for 2D materials
    
    # Add legend if we plotted anything
    if plotted_surfaces > 0 or len(iso_values) > 0:
        # Create custom legend
        legend_elements = []
        for iso in iso_values:
            if cube.data_type == 'density' and iso < 0:
                continue
            color = colors_map(norm(iso))
            from matplotlib.lines import Line2D
            legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor=color, markersize=10,
                                        label=f'iso = {iso:.2e}'))
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))
    
    plt.tight_layout()
    return fig, ax


def plot_slice(cube: CubeFile, axis: str = 'z', position: Optional[int] = None,
               cmap: Optional[str] = None, show_contours: bool = True,
               fig_ax: Optional[Tuple] = None) -> Tuple:
    """Plot 2D slice of cube data."""
    
    if fig_ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
    else:
        fig, ax = fig_ax
    
    # Get slice data
    slice_data = cube.get_slice(axis, position)
    
    # Choose colormap
    if cmap is None:
        cmap = {
            'density': 'viridis',
            'spin': 'RdBu',
            'potential': 'RdBu'
        }.get(cube.data_type, 'viridis')
    
    # Determine extent
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
    if cube.data_type == 'spin' or cube.data_type == 'potential':
        vmax = np.max(np.abs(slice_data))
        vmin = -vmax
    else:
        vmin, vmax = np.min(slice_data), np.max(slice_data)
    
    im = ax.imshow(slice_data.T, origin='lower', extent=extent, cmap=cmap,
                   vmin=vmin, vmax=vmax, aspect='auto')
    
    # Add contours
    if show_contours and np.std(slice_data) > 1e-10:
        X = np.linspace(extent[0], extent[1], slice_data.shape[0])
        Y = np.linspace(extent[2], extent[3], slice_data.shape[1])
        X, Y = np.meshgrid(X, Y)
        
        if cube.data_type == 'density':
            levels = np.logspace(np.log10(np.max(slice_data) * 0.01), 
                               np.log10(np.max(slice_data) * 0.9), 8)
        else:
            levels = 10
        
        contours = ax.contour(X, Y, slice_data.T, levels=levels, colors='black', 
                             alpha=0.3, linewidths=0.5)
        ax.clabel(contours, inline=True, fontsize=8, fmt='%.2e')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar_label = {
        'density': 'Density (e/Bohr³)',
        'spin': 'Spin density (e/Bohr³)', 
        'potential': 'Potential (Hartree)'
    }.get(cube.data_type, 'Value')
    cbar.set_label(cbar_label, fontsize=12)
    
    # Set labels and title
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    
    pos_value = position if position is not None else cube.nvoxels[['x','y','z'].index(axis)] // 2
    title = f"{cube.data_type.capitalize()} - Slice along {axis.upper()} at index {pos_value}"
    ax.set_title(title, fontsize=14)
    
    plt.tight_layout()
    return fig, ax


def interactive_menu(cube: CubeFile):
    """Interactive menu for choosing visualization options."""
    print("\n" + "="*60)
    print("VISUALIZATION OPTIONS")
    print("="*60)
    
    # Show recommended isovalues
    print(f"\nRecommended isovalues for {cube.data_type}:")
    print(f"  {', '.join([f'{x:.2e}' for x in cube.recommended_iso])}")
    
    # Ask user for choice
    print("\nChoose visualization:")
    print("1. 3D isosurfaces with recommended values")
    print("2. 3D isosurfaces with custom values")
    print("3. 2D slice through the data")
    print("4. Multi-panel view (3D + slices)")
    print("5. Just print data statistics")
    
    choice = input("\nEnter choice (1-5) [1]: ").strip() or "1"
    
    if choice == "1":
        return "3d", cube.recommended_iso, None
    
    elif choice == "2":
        print("\nEnter custom isovalues (comma-separated):")
        print("Example: -0.01,-0.001,0.001,0.01")
        iso_input = input("Values: ").strip()
        try:
            iso_values = [float(x.strip()) for x in iso_input.split(',')]
            return "3d", iso_values, None
        except:
            print("Invalid input. Using recommended values.")
            return "3d", cube.recommended_iso, None
    
    elif choice == "3":
        axis = input("Slice axis (x/y/z) [z]: ").strip().lower() or "z"
        try:
            position = input(f"Slice position (0-{cube.nvoxels[['x','y','z'].index(axis)]-1}) [middle]: ").strip()
            position = int(position) if position else None
        except:
            position = None
        return "slice", None, (axis, position)
    
    elif choice == "4":
        return "multi", cube.recommended_iso, None
    
    else:
        return "stats", None, None


def create_multi_panel_figure(cube: CubeFile, iso_values: Optional[List[float]] = None):
    """Create a figure with multiple panels showing different views."""
    fig = plt.figure(figsize=(20, 12))
    
    # 3D isosurface
    ax1 = fig.add_subplot(221, projection='3d')
    plot_isosurface(cube, iso_values, fig_ax=(fig, ax1))
    
    # XY slice (through middle of Z)
    ax2 = fig.add_subplot(222)
    plot_slice(cube, axis='z', fig_ax=(fig, ax2))
    
    # XZ slice (through middle of Y)
    ax3 = fig.add_subplot(223)
    plot_slice(cube, axis='y', fig_ax=(fig, ax3))
    
    # Data statistics
    ax4 = fig.add_subplot(224)
    ax4.axis('off')
    
    # Create text summary
    stats_text = f"""Data Statistics:
    
File: {os.path.basename(cube.filename)}
Type: {cube.data_type}
Structure: {'2D material' if cube.is_2d else '3D material'}

Grid: {cube.nvoxels[0]} × {cube.nvoxels[1]} × {cube.nvoxels[2]}
Total points: {np.prod(cube.nvoxels):,}

Data range: [{np.min(cube.data):.6f}, {np.max(cube.data):.6f}]
Mean: {np.mean(cube.data):.6f}
Std: {np.std(cube.data):.6f}

Non-zero points: {np.sum(cube.data != 0):,} ({100*np.sum(cube.data != 0)/cube.data.size:.1f}%)
"""
    
    if cube.data_type == 'spin':
        total_up = np.sum(cube.data[cube.data > 0])
        total_down = abs(np.sum(cube.data[cube.data < 0]))
        stats_text += f"\nSpin-up excess: {total_up:.6f}"
        stats_text += f"\nSpin-down excess: {total_down:.6f}"
        stats_text += f"\nNet spin: {np.sum(cube.data):.6f}"
    
    ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, 
             fontsize=12, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    return fig


def main():
    """Main function to handle command line arguments and plotting."""
    parser = argparse.ArgumentParser(description='Visualize cube files from CRYSTAL calculations')
    parser.add_argument('filename', help='Cube file to visualize')
    parser.add_argument('--iso', type=str, help='Override isosurface values (comma-separated)')
    parser.add_argument('--cmap', help='Colormap name')
    parser.add_argument('--alpha', type=float, default=0.7, help='Transparency (0-1)')
    parser.add_argument('--slice', nargs=2, help='Slice axis and position')
    parser.add_argument('--no-atoms', action='store_true', help='Do not show atoms')
    parser.add_argument('--save', help='Save figure to file')
    parser.add_argument('--view', type=str, help='View angles as elev,azim')
    parser.add_argument('--interactive', '-i', action='store_true', 
                      help='Interactive mode with menu')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.filename):
        print(f"Error: File '{args.filename}' not found!")
        sys.exit(1)
    
    # Read cube file
    print(f"Reading cube file: {args.filename}")
    cube = CubeFile(args.filename)
    cube.print_info()
    
    # Handle all-zero data
    if np.all(cube.data == 0):
        print("\nWARNING: This file contains all zeros!")
        print("This is common for:")
        print("- Non-magnetic materials (SPIN files)")
        print("- Systems with no charge transfer")
        if not args.interactive:
            print("\nNothing to visualize. Exiting.")
            sys.exit(0)
    
    # Interactive mode
    if args.interactive:
        viz_type, iso_values, slice_params = interactive_menu(cube)
        
        if viz_type == "stats":
            # Already printed
            sys.exit(0)
        elif viz_type == "3d":
            fig, ax = plot_isosurface(cube, iso_values, cmap=args.cmap, 
                                    alpha=args.alpha, show_atoms=not args.no_atoms)
        elif viz_type == "slice":
            fig, ax = plot_slice(cube, axis=slice_params[0], position=slice_params[1], 
                               cmap=args.cmap)
        elif viz_type == "multi":
            fig = create_multi_panel_figure(cube, iso_values)
    
    # Command line mode
    else:
        # Parse isosurface values
        iso_values = None
        if args.iso:
            try:
                iso_values = [float(x.strip()) for x in args.iso.split(',')]
            except:
                print("Warning: Invalid iso values. Using recommended values.")
        
        # Parse view angles
        view_angles = None
        if args.view:
            try:
                elev, azim = map(float, args.view.split(','))
                view_angles = (elev, azim)
            except:
                print("Warning: Invalid view angles. Using default.")
        
        # Create plots
        if args.slice:
            # 2D slice plot
            axis = args.slice[0].lower()
            try:
                position = int(args.slice[1])
            except:
                position = None
            fig, ax = plot_slice(cube, axis=axis, position=position, cmap=args.cmap)
        else:
            # 3D isosurface plot
            fig, ax = plot_isosurface(cube, iso_values=iso_values, cmap=args.cmap,
                                    alpha=args.alpha, show_atoms=not args.no_atoms,
                                    view_angles=view_angles)
    
    # Save or show figure
    if args.save:
        plt.savefig(args.save, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
