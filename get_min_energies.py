import glob, os, operator
from matplotlib import pyplot as plt

# Gets the total energy from each output file
def get_energy(file):
    with open(file,'r') as f:
        for i,line in enumerate(f):
            if 'TOTAL ENERGY + DISP' in line:
                energy = float(line.split()[-1])
    return(energy)

# Goes through all files in directory, grabs energies and sorts them
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

sorted_d = sorted(dict(zip(materials,energies)).items(), key=operator.itemgetter(1))
first3pairs = sorted_d[:3] #print 3 min energies
print(first3pairs)

# Quick plotting
plt.figure(figsize=(10,5),dpi=300)
for pair in sorted_d:
    plt.scatter(pair[0],pair[1])
plt.xticks(rotation=89)
plt.ylabel('Energy (Hartree)')
plt.tight_layout()
plt.savefig('energy_plot.png',format='png')
#plt.show()
