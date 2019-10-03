#
#
import sys
from Bio.PDB.PDBParser import PDBParser
import matplotlib.pyplot as plt
import numpy as np
from scipy import  stats as stats

parser = PDBParser(PERMISSIVE=1, QUIET=True)
name = "6AAC"
folder = "PED6AAC-k18_pool"

Rg_array_pool = []
infile = open(folder[3:]+"_Rg_values.txt", 'r')
for line in infile:
	line = line.strip()
	if len(line)>1:
		Rg = float(line)
		Rg_array_pool.append(Rg)
infile.close()
'''
for i in range(0, 20000):
	filename = folder+"/"+str(i)+"a_130.pdb"
	structure = parser.get_structure(name+"-"+str(i), filename)
	atoms = structure[0].get_atoms()
	x = 0.0
	y = 0.0
	z = 0.0
	N = 0.0
	for atom in atoms:
		x += atom.get_coord()[0]
		y += atom.get_coord()[1]
		z += atom.get_coord()[2]
		N += 1
	centr = []
	centr.append( x/N )
	centr.append( y/N )
	centr.append( z/N )
	sumD = 0.0
	for atom2 in structure[0].get_atoms():
		x = atom2.get_coord()[0]
		y = atom2.get_coord()[1]
		z = atom2.get_coord()[2]
		sumD += np.sqrt( np.square(x-centr[0]) + np.square(y-centr[1]) + np.square(z-centr[2]) )
	Rg = sumD / N
	print Rg
	Rg_array_pool.append( Rg )
'''
#------
Rg_array = []
infile = open(name+"_Rg_values.txt", 'r')
for line in infile:
	line = line.strip()
	if len(line)>1:
		Rg = float(line)
		Rg_array.append(Rg)

infile.close()


t, p = stats.ttest_ind(Rg_array, Rg_array_pool)
p = "%.2E" % p
p1 = plt.hist(Rg_array, bins=9, histtype='stepfilled', alpha=0.5, normed=True, color='r')
p2 = plt.hist(Rg_array_pool, bins=9, histtype='stepfilled', alpha=0.35, normed=True, color='b')
plt.suptitle("Rg distribution of "+name+" ensemble vs. its pool")
plt.title("P(ensemble,pool)="+p, fontsize=10)
plt.xlabel("Rg (Angstrom)")
plt.ylabel("density")
plt.ylim(0, 0.08)
plt.xlim(0,70)
plt.legend(("ensemble", "pool"), fontsize=12)
plt.savefig("Rg_"+name+"_vs_pool.png", dpi=600)
