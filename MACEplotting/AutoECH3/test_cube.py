#!/usr/bin/env python3
"""
Quick test script to check if a cube file loads correctly.
"""

import sys
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python test_cube.py filename.cube")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, 'r') as f:
        # Read header
        comment1 = f.readline().strip()
        comment2 = f.readline().strip()
        
        # Read atoms and origin
        line = f.readline().split()
        natoms = int(line[0])
        origin = [float(line[1]), float(line[2]), float(line[3])]
        
        # Read voxel information
        nvoxels = []
        for i in range(3):
            line = f.readline().split()
            nvoxels.append(int(line[0]))
        
        print(f"File: {filename}")
        print(f"Comments: {comment1} | {comment2}")
        print(f"Number of atoms: {abs(natoms)}")
        print(f"Grid size: {nvoxels[0]} x {nvoxels[1]} x {nvoxels[2]}")
        print(f"Total grid points: {nvoxels[0] * nvoxels[1] * nvoxels[2]}")
        
        # Count data points to read a sample
        total_points = nvoxels[0] * nvoxels[1] * nvoxels[2]
        
        # Skip atoms
        for i in range(abs(natoms)):
            f.readline()
        
        # Read first few hundred data points to get statistics
        sample_data = []
        lines_to_read = min(100, total_points // 10)  # Read ~10% or 100 lines
        for i in range(lines_to_read):
            line = f.readline()
            if line:
                sample_data.extend([float(x) for x in line.split()])
        
        if sample_data:
            sample_data = np.array(sample_data)
            print(f"\nData statistics (from sample):")
            print(f"  Min value: {np.min(sample_data):.6f}")
            print(f"  Max value: {np.max(sample_data):.6f}")
            print(f"  Mean: {np.mean(sample_data):.6f}")
            print(f"  Std dev: {np.std(sample_data):.6f}")
        
        # Quick check of data type from filename
        if 'spin' in filename.lower():
            print("\nData type: Spin density")
            if sample_data is not None:
                max_abs = np.max(np.abs(sample_data))
                if max_abs < 0.01:
                    print("Suggested isovalues: -0.001, -0.0001, 0.0001, 0.001")
                elif max_abs < 0.1:
                    print("Suggested isovalues: -0.01, -0.001, 0.001, 0.01")
                else:
                    print("Suggested isovalues: -0.05, -0.01, 0.01, 0.05")
            else:
                print("Suggested isovalues: -0.01, -0.001, 0.001, 0.01")
        elif 'ech3' in filename.lower() or ('density' in filename.lower() and 'spin' not in filename.lower()):
            print("\nData type: Likely electron density")
            print("Suggested isovalues: 0.001, 0.01, 0.05")
        elif 'pot' in filename.lower() or 'potential' in filename.lower():
            print("\nData type: Likely electrostatic potential")
            print("Suggested isovalues: -0.05, -0.01, 0.01, 0.05")
        
        print("\nFile appears to be valid!")
        
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)
