import numpy as np
import os, glob, shutil
from collections import Counter, OrderedDict

# Get a dictionary with the number of atoms and atom types for a cif file
def numatoms(file):
    with open(file, 'r') as f:
        count = 9999999
        list = []
        for i, line in enumerate(f):
            if line.startswith('_atom_site_fract_z'):
                count=i
            if i >= count+1:
                list.append(line.split()[0])
        return(dict(Counter(list)))

# Group materials with the same stoichiometry and make list
def group_unique_materials(list_dict,unique_compounds):
    groups = [[] for x in range(len(unique_compounds))]
    for i in range(len(unique_compounds)):
        for j in list_dict:
            if j[1] == unique_compounds[i]:
                groups[i].append(j[0])
    #print(groups)
    return(groups)

# Make directories and copy files with the same stoichiometry
def make_dirs_copy_files(DIR,unique_compounds,group_list):
    for i in range(len(unique_compounds)):
        os.mkdir(os.path.join(DIR+str(unique_compounds[i])))
        for j in group_list[i]:
            shutil.copy(os.path.join(DIR+str(j)+'.cif'),os.path.join(DIR+str(unique_compounds[i])))
            shutil.copy(os.path.join(DIR+str(j)+'.out'),os.path.join(DIR+str(unique_compounds[i])))
            shutil.copy(os.path.join(DIR+str(j)+'.d12'),os.path.join(DIR+str(unique_compounds[i])))


# Parses through all CIF files in directory // Gets unique compounds and groups them in a nested list
DIR = (os.getcwd()+'/')
nDIR = len(DIR)
ntype = len(".cif")
pathlist = glob.glob(DIR+'*.cif')
list_dict = []
compositions = []
for path in pathlist:
    path_in_str = str(path)
    material = path_in_str[nDIR:-ntype]
    # Make a list like this: [[material1, {ordered elements}], [material2, {ordered elements}], ...]
    # Ordered dictionaries can be compared directly, so we can group materials with same stoichiometry
    list_dict.append([material,str(OrderedDict(sorted(numatoms(path_in_str).items())))])
    compositions.append(str(OrderedDict(sorted(numatoms(path_in_str).items()))))
unique_compounds = list(set(compositions))
group_list = group_unique_materials(list_dict,unique_compounds)
make_dirs_copy_files(DIR,unique_compounds,group_list)
#print(unique_compounds)
#print(group_list)