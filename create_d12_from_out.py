 
#!/usr/bin/env python

"""
Takes optimized geometry from CRYSTAL17/CRYSTAL23 output and edits the input .d12 file with the new geometry.
Usage: create_d12_from_optout.py
Danny Maldonado-Lopez
"""

import os, sys, math, subprocess
import numpy as np

subprocess.run('sed -i "/^ PROCESS/d" *.out',shell=True, check=True)


def get_geom(material):
    spacegroup = 9999
    cryin = open(material+".d12")
    for i, line in enumerate(cryin):
        if i==1 and line.startswith('CRYSTAL'):
            struc = 'bulk'
        elif i==1 and line.startswith('SLAB'):
            struc = 'slab'
        if i==2 and struc == 'slab':
            spacegroup = int(line)
        if i==3 and struc == 'bulk':
            spacegroup = int(line)
        if line.startswith('SLABCUT'):
            struc = 'slabcut'
    cryin.close()

    cryout = open(material+".out")
    md = 999999
    atomnum = np.zeros(999)
    atomx = np.zeros(999)
    atomy = np.zeros(999)
    atomz = np.zeros(999)
    unique = []
    print(struc)

    if struc == "bulk":
        for i, line in enumerate(cryout):
            if "CRYSTAL FAMILY" in line:
                crystal_fam = str(line.split()[-1])
            if line.startswith(" NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL:"):
                natoms_sym = int(line.split()[-1])
            if "FINAL OPTIMIZED GEOMETRY" in line:
                md = i
                continue
            if i == md+6:
                cell    = line.split()
                a = float(cell[0])
                b = float(cell[1])
                c = float(cell[2])
                alpha = float(cell[3])
                beta = float(cell[4])
                gamma = float(cell[5])
            if i == md+8:
                natoms = int(line.split()[-1])
            if md+11 <= i and i<= md+10+natoms:
                atomnum[i-md-11] = int(line.split()[-5])
                unique.append(str(line.split()[1]))
                atomx[i-md-11] = float(line.split()[-3])
                atomy[i-md-11] = float(line.split()[-2])
                atomz[i-md-11] = float(line.split()[-1])
        unique_atoms=[]
        for i in range(0,natoms):
            if unique[i] == "T":
                unique_atoms.append([str(int(atomnum[i]))+'   '+str(atomx[i])+'   '+str(atomy[i])+'   '+str(atomz[i])])
        # Lattice parameters to be printed depend on crystal family
        if crystal_fam == "TRICLINIC":
            latpar = str(a)+'   '+str(b)+'   '+str(c)+'   '+str(alpha)+'   '+str(beta)+'   '+str(gamma)
        if crystal_fam == "MONOCLINIC":
            latpar = str(a)+'   '+str(b)+'   '+str(c)+'   '+str(beta)
        if crystal_fam == "ORTHORHOMBIC":
            latpar = str(a)+'   '+str(b)+'   '+str(c)
        if crystal_fam == "TETRAGONAL":
            latpar = str(a)+'   '+str(c)
        if crystal_fam == "TRIGONAL":
            if a == b == c:
                latpar = str(a)+'   '+str(alpha)
            else:
                latpar = str(a)+'   '+str(c)
        if crystal_fam == "HEXAGONAL":
            latpar = str(a)+'   '+str(c)
        if crystal_fam == "CUBIC":
            latpar = str(a)

    if struc == "slab":
        for i, line in enumerate(cryout):
            if line.startswith(" TWO-SIDED PLANE GROUP"):
                spacegroup = int(line.split()[4])
            if line.startswith(" NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL:"):
                natoms_sym = int(line.split()[-1])
            if "FINAL OPTIMIZED GEOMETRY" in line:
                md = i
                continue
            if i == md+6:
                cell    = line.split()
                a = float(cell[0])
                b = float(cell[1])
                c = float(cell[2])
                alpha = float(cell[3])
                beta = float(cell[4])
                gamma = float(cell[5])
            if i == md+8:
                natoms = int(line.split()[-1])
            if md+11 <= i and i<= md+10+natoms:
                atomnum[i-md-11] = int(line.split()[-5])
                unique.append(str(line.split()[1]))
                atomx[i-md-11] = float(line.split()[-3])
                atomy[i-md-11] = float(line.split()[-2])
                atomz[i-md-11] = float(line.split()[-1])
        unique_atoms=[]
        for i in range(0,natoms):
            if unique[i] == "T":
                unique_atoms.append([str(int(atomnum[i]))+'   '+str(atomx[i])+'   '+str(atomy[i])+'   '+str(atomz[i])])
        # For slabs, we can just all lattice paramters, without issues
        latpar=str(a)+'   '+str(b)+'   '+str(gamma)

    # For slabcut, I am exporting without symmetry
    if struc == "slabcut":
        for i, line in enumerate(cryout):
            if "FINAL OPTIMIZED GEOMETRY" in line:
                md = i
                continue
            if i == md+6:
                cell    = line.split()
                a = float(cell[0])
                b = float(cell[1])
                c = float(cell[2])
                alpha = float(cell[3])
                beta = float(cell[4])
                gamma = float(cell[5])
            if i == md+8:
                natoms = int(line.split()[-1])
                natoms_sym = natoms
            if md+11 <= i and i<= md+10+natoms:
                atomnum[i-md-11] = int(line.split()[-5])
                unique.append(str(line.split()[1]))
                atomx[i-md-11] = float(line.split()[-3])
                atomy[i-md-11] = float(line.split()[-2])
                atomz[i-md-11] = float(line.split()[-1])
        unique_atoms=[]
        for i in range(0,natoms):
            unique_atoms.append([str(int(atomnum[i]))+'   '+str(atomx[i])+'   '+str(atomy[i])+'   '+str(atomz[i])])
        latpar=str(a)+'   '+str(b)+'   '+str(gamma)

    return(struc,latpar,natoms_sym,unique_atoms)

def print_d12(material,latpar,natoms_sym,unique_atoms):
    cryin=open(material+".d12")
    if struc == 'bulk':
        endcounter = 0
        lines = []
        for i,line in enumerate(cryin):
            lines.append(line)
            if line.startswith('OPTGEOM'):
                optgeom_line_number = i
            if line.startswith('END'):
                endcounter += 1
            if endcounter == 1:
                end_line_number = i+1
        lines[4]=str(latpar)+'\n'
        lines[5]=str(natoms_sym)+'\n'
        linecounter = 6
        for i in range(len(unique_atoms)):
            lines[linecounter+i]=unique_atoms[i][0]+'\n'
        lines = lines[:max(optgeom_line_number, 0)] + lines[end_line_number:]
    cryin.close()

    cryin=open(material+".d12")
    if struc == 'slab':
        endcounter = 0
        lines = []
        for i,line in enumerate(cryin):
            lines.append(line)
            if line.startswith('OPTGEOM'):
                optgeom_line_number = i
            if line.startswith('END'):
                endcounter += 1
            if endcounter == 1:
                end_line_number = i+1
        lines[3]=str(latpar)+'\n'
        lines[4]=str(natoms_sym)+'\n'
        linecounter = 5
        for i in range(len(unique_atoms)):
            lines[linecounter+i]=unique_atoms[i][0]+'\n'
        lines = lines[:max(optgeom_line_number, 0)] + lines[end_line_number:]
    cryin.close()

    cryin=open(material+".d12")
    if struc == 'slabcut':
        endcounter = 0
        lines = []
        for i,line in enumerate(cryin):
            lines.append(line)
            if line.startswith('SLABCUT'):
                slabcut_line_number = i
            if line.startswith('END'):
                endcounter += 1
            if endcounter == 1:
                end_line_number = i+1
        lines[1]='SLAB\n'
        lines[2]='1\n'
        lines[3]=str(latpar)+'\n'
        lines[4]=str(natoms_sym)+'\n'
        linecounter = 5
        lines = lines[:max(linecounter+1, 0)] + lines[end_line_number-1:]
        for i in range(len(unique_atoms)):
            lines[linecounter+i]=unique_atoms[i][0]+'\n'
    cryin.close()

    new_cryin = open(material+'_sp.d12','w+')
    for item in lines:
        new_cryin.write(item)
    new_cryin.close()


data_files = os.listdir(os.getcwd())
for file_name in data_files:
    if ".d12" in file_name:
        file_name_str = str(file_name)
        material = file_name_str[:-4]
        print(material)
        struc,latpar,natoms_sym,unique_atoms=get_geom(material)
        print_d12(material,latpar,natoms_sym,unique_atoms)
