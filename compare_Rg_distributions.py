#
#
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import sys

def calc_Rg_distr(filename, Npairs_res, name):
	myhash = {}
	infile = open(filename, 'r')
	for line in infile:
		line = line.strip()
		if len(line)>1:
			array = line.split(";")
			pair = ( int(array[0]), int(array[1]) )
			numbers = array[2][1:-1]
			number_array = numbers.split(", ")
			num_array = []
			for i in number_array:
				num_array.append(float(i))
			if pair in myhash.keys():
				myhash[pair] += num_array
			else:
				myhash[pair] = num_array
	infile.close()

	n_conf = len(myhash[(1,2)])
	Npairs_conf = (n_conf*(n_conf-1))/2
	Rgs = []
	for c in range(0,n_conf):
		dist_jk = 0
		for j in range(1,last_res ):
			for k in range(j+1, last_res+1):
				dist_jk += ( myhash[(j,k)][c]**2 )
		Rg = np.sqrt( dist_jk / Npairs_res )
		Rgs.append(Rg)
	print name+":  <Rg>="+str(np.mean(Rgs))+" A"
	return Rgs


ens_type = "ensemble" # "ensemble" or "random_pool"
res_folder = "results/"+ens_type+"/"
ensemble = "PED00001e001"
ens_type2 = "ensemble" # "ensemble" or "random_pool"
ensemble2 = "PED00001e002"
atomtype = "CA"
first_res = -1
last_res = 90

Npairs_res = (last_res*(last_res-1))/2
filename1 = "results/"+ens_type+"/"+ensemble+"_"+atomtype+"_distance_distributions.txt"
filename2 = "results/"+ens_type2+"/"+ensemble2+"_"+atomtype+"_distance_distributions.txt"
Rg_array1 = calc_Rg_distr(filename1, Npairs_res, ensemble)
print "vs"
Rg_array2 = calc_Rg_distr(filename2, Npairs_res, ensemble2)

p = stats.mannwhitneyu(Rg_array1, Rg_array2)[1]
p = "%.2E" % p
p1 = plt.hist(Rg_array1, bins=9, histtype='stepfilled', alpha=0.5, normed=True, color='r')
p2 = plt.hist(Rg_array2, bins=9, histtype='stepfilled', alpha=0.35, normed=True, color='b')
plt.suptitle("Rg distributions of "+ensemble+" & "+ensemble2)
plt.title("Difference(U-test):  p="+p, fontsize=10)
plt.xlabel("Rg (Angstrom)")
plt.ylabel("density")
plt.legend((ensemble, ensemble2), fontsize=12)
plt.savefig(res_folder+"Rg_"+ensemble+"-"+ensemble2+".png", dpi=600)
