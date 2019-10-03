import Bio.PDB
import numpy
import sys

window_sizes = ['3','5','7']

pdb_code = "PED6AAC-k18_pool"
#pdb_code = "PED7AAC-ntail_pool"
pdb_filename = "PED6AAC-k18_pool/ensemble_small.pdb"
#pdb_filename = "PED7AAC-ntail_pool/ensemble_small.pdb"
super_imposer = Bio.PDB.Superimposer()
structure = Bio.PDB.PDBParser(QUIET=True).get_structure(pdb_code, pdb_filename)
prot_length = 130
if "7AAC" in pdb_code:
	prot_length = 132

for size in window_sizes:
	outfile = open(pdb_code+"_local_rmsd_"+size+"aa-S.txt", 'w')
	size = int(size)
	for i in range(1, prot_length+1-size ):
		print "Calculating superposition for aa. "+str(i)+"-"+str(i+size-1)+" of "+pdb_code+" with "+str(size)+"aa window"
		rmsd_array = []
		ref_model = structure[0]
		ref_atoms = []
		for ref_chain in ref_model:
			for ref_res in ref_chain:
				if ref_res.get_id()[1] >= i and ref_res.get_id()[1] < i+size:
					ref_atoms.append(ref_res['CA'])
					ref_atoms.append(ref_res['N'])
					ref_atoms.append(ref_res['C'])
					ref_atoms.append(ref_res['O'])

		for j in range(1,len(structure)):
			#Build paired lists of c-alpha atoms, ref_atoms and alt_atoms
			alt_atoms = []
			for alt_chain in structure[j]:
				for alt_res in alt_chain :
					if alt_res.get_id()[1] >= i and alt_res.get_id()[1] < i+size:
						alt_atoms.append(alt_res['CA'])              
						alt_atoms.append(alt_res['N'])               
						alt_atoms.append(alt_res['C'])            
						alt_atoms.append(alt_res['O'])

			#Align these paired atom lists:
			super_imposer.set_atoms(ref_atoms, alt_atoms)
			rmsd_array.append( super_imposer.rms )

		outfile.write(str(i)+";"+str(i+size-1)+";"+str(rmsd_array)+"\n")
	print

	outfile.close()


