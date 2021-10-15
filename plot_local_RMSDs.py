#
#
import numpy as np
import matplotlib.pyplot as plt
import sys

ens_type = "ensemble" # "ensemble" or "random_pool"
res_folder = "results/"+ens_type+"/"
ensemble_file = "PED00001e001.pdb"
first_res = -1
last_res = 90


infile = open(res_folder+ensemble_file[0:-4]+"_local_rmsd_5aa.txt", 'r')
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

x_pos = np.arange(first_res, last_res+1+1-5)
ax = plt.subplot(1,1,1)
plt.fill_between(x_pos, uppers, lowers, color='gray', alpha=0.6)
plt.plot(x_pos, medians, 'k', label="medians")
plt.plot(x_pos, uppers, 'gray', label="95% percentiles")
plt.xlim(first_res, last_res)
plt.ylim(0, 6)
plt.title("Local superimposability of "+ensemble_file[0:-4])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)
plt.xlabel("residue number")
plt.ylabel("5-residue local backbone RMSD (Angstrom)")
plt.savefig(res_folder+ensemble_file[0:-4]+"_local_rmsd_5aa.png", dpi=600)
plt.clf()
