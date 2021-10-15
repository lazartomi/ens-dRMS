import Bio.PDB

ens_type = "ensemble" # "ensemble" or "random_pool"
data_folder = "data/"+ens_type+"/"
res_folder = "results/"+ens_type+"/"
ensemble_file = "PED00001e001.pdb"
first_res = -1
last_res = 90

super_imposer = Bio.PDB.Superimposer()
structure = Bio.PDB.PDBParser(QUIET=True).get_structure(ensemble_file[0:-4], data_folder+ensemble_file)
prot_length = last_res - first_res

outfile = open(res_folder+ensemble_file[0:-4]+"_local_rmsd_5aa.txt", 'w')
for i in range(first_res, prot_length+1-5 ):
	#print "Calculating superposition for aa. "+str(i)+"-"+str(i+5-1)+" of "+ensemble_file[0:-4]+" with 5aa window"
	rmsd_array = []
	ref_model = structure[0]
	ref_atoms = []
	for ref_chain in ref_model:
		for ref_res in ref_chain:
			if ref_res.get_id()[1] >= i and ref_res.get_id()[1] < i+5:
				ref_atoms.append(ref_res['CA'])
				ref_atoms.append(ref_res['N'])
				ref_atoms.append(ref_res['C'])
				if 'O' in ref_res:
					ref_atoms.append(ref_res['O'])

	for j in range(1,len(structure)):
		#Build paired lists of c-alpha atoms, ref_atoms and alt_atoms
		alt_atoms = []
		for alt_chain in structure[j]:
			for alt_res in alt_chain :
				if alt_res.get_id()[1] >= i and alt_res.get_id()[1] < i+5:
					alt_atoms.append(alt_res['CA'])              
					alt_atoms.append(alt_res['N'])               
					alt_atoms.append(alt_res['C']) 
					if 'O' in alt_res:
						alt_atoms.append(alt_res['O'])

		#Align these paired atom lists:
		super_imposer.set_atoms(ref_atoms, alt_atoms)
		rmsd_array.append( super_imposer.rms )

	outfile.write(str(i)+";"+str(i+5-1)+";"+str(rmsd_array)+"\n")

outfile.close()
