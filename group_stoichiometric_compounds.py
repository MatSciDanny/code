import os, glob, shutil, operator
from collections import Counter, OrderedDict
from matplotlib import pyplot as plt

# Get a dictionary with the number of atoms and atom types from a cif file
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
    return(groups)

# Make directories and copy files with the same stoichiometry
def make_dirs_copy_files(DIR,unique_compounds,group_list):
    all_dirs = []
    for i in range(len(unique_compounds)):
        dir_name = str(unique_compounds[i]).replace("[","").replace("]","").replace("(","").replace(")","").replace("'","").replace("OrderedDict","").replace(",","_").replace(" ","")
        all_dirs.append(dir_name)
        os.mkdir(os.path.join(DIR+str(dir_name)))
        for j in group_list[i]:
            shutil.copy(os.path.join(DIR+str(j)+'.cif'),os.path.join(DIR+dir_name))
            shutil.copy(os.path.join(DIR+str(j)+'.out'),os.path.join(DIR+dir_name))
            shutil.copy(os.path.join(DIR+str(j)+'.d12'),os.path.join(DIR+dir_name))
    return(all_dirs)

# Get total energy from each output file
def get_energy(file):
    with open(file,'r') as f:
        for i,line in enumerate(f):
            if 'TOTAL ENERGY + DISP' in line:
                energy = float(line.split()[-1])
    return(energy)

# Go into folders and plot energy vs file name (from low to high energy) for each stoichiometry
def sort_plot_energies(all_dirs):
    for i in all_dirs:
        parent_dir = os.getcwd()
        os.chdir(str(i))
        DIR = (os.getcwd()+'/')
        nDIR = len(DIR)
        ntype = len(".out")
        pathlist = glob.glob(DIR+'*.out')
        materials,energies = [],[]
        for path in pathlist:
            path_in_str = str(path)
            print(path_in_str)
            material = path_in_str[nDIR:-ntype]
            materials.append(material)
            energies.append(get_energy(path_in_str))
        #Sort directories from low to high energy
        sorted_d = sorted(dict(zip(materials,energies)).items(), key=operator.itemgetter(1))
        first3pairs = sorted_d[:3] #print 3 min energies
        print(first3pairs)
        # Plotting
        plt.figure(figsize=(10,5),dpi=300)
        for pair in sorted_d:
            plt.scatter(pair[0],pair[1])
        plt.xticks(rotation=89)
        plt.ylabel('Energy (Hartree)')
        plt.tight_layout()
        plt.savefig('energy_plot.png',format='png')
        os.chdir(parent_dir)


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
    # Ordered dictionaries can be compared directly to group materials with same stoichiometry
    list_dict.append([material,str(OrderedDict(sorted(numatoms(path_in_str).items())))])
    compositions.append(str(OrderedDict(sorted(numatoms(path_in_str).items()))))
unique_compounds = list(set(compositions))
# Now we use the functions we defined previously
group_list = group_unique_materials(list_dict,unique_compounds)
all_dirs = make_dirs_copy_files(DIR,unique_compounds,group_list)
sort_plot_energies(all_dirs)