#
#
import numpy as np
import matplotlib.pyplot as plt
import sys

window_sizes = ['3','5','7']

pdb_code = "PED6AAC-k18_pool"
#pdb_code = "PED7AAC-ntail_pool"

for window_size in window_sizes:
	infile = open(pdb_code+"_local_rmsd_"+window_size+"aa-S.txt", 'r')
	medians = []
	uppers = []
	lowers = []
	stdevs = []
	for line in infile:
		line = line.strip()
		if len(line)>1:
			array = line.split(";")
			rmsd_array = array[2][1:-1].split(", ")
			array2 = []
			for rmsd in rmsd_array:
				array2.append( float(rmsd) )
			rmsd_array = array2
			medians.append( np.median(rmsd_array) )
			stdevs.append( np.std(rmsd_array) )
			uppers.append( np.percentile(rmsd_array, 95) )
			lowers.append( np.percentile(rmsd_array, 5) )
	infile.close()

	x_pos = np.arange(len(medians))
	ax = plt.subplot(1,1,1)
	plt.fill_between(x_pos, uppers, lowers, color='gray', alpha=0.6)
	plt.plot(x_pos, medians, 'k', label="medians")
	plt.plot(x_pos, uppers, 'gray', label="95% percentiles")
	plt.xlim(1,len(medians)-1)
	plt.ylim(0, 6)
	plt.title("Local superimposability of "+pdb_code)
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles, labels)
	plt.xlabel("residue number")
	plt.ylabel(window_size+"-residue local backbone RMSD (Angstrom)")
	plt.savefig(pdb_code+"_local_rmsd_"+window_size+"aa-S.png", dpi=600)
	plt.clf()
