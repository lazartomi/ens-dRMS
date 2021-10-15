#
#
import scipy.stats as stat
import numpy as np

def read_distr(filename):
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
			myhash[pair] = num_array
	infile.close()
	return myhash



ens_type = "ensemble" # "ensemble" or "random_pool"
res_folder = "results/"+ens_type+"/"
ensemble = "PED00001e001"
ensemble2 = "PED00001e002"
atomtype = "CA"
first_res = -1
last_res = 90


filename1 = "results/"+ens_type+"/"+ensemble+"_"+atomtype+"_distance_distributions.txt"
myhash = read_distr(filename1)
filename2 = "results/"+ens_type+"/"+ensemble2+"_"+atomtype+"_distance_distributions.txt"
myhash2 = read_distr(filename2)

outfilename = "results/"+ens_type+"/"+ensemble+"-"+ensemble2+"_"+atomtype+"_dist_distrib_differences.txt"
outfile = open(outfilename, 'w')
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			outfile.write(str(j)+";"+str(k)+";1.0;0.0;0.0;0.0\n")
		else:
			pair = (j,k)
			median_diff = abs( np.median(myhash[pair])-np.median(myhash2[pair]) )
			std_diff = abs( np.std(myhash[pair])-np.std(myhash2[pair]) )
			bigger = np.median(myhash[pair])
			if np.median(myhash2[pair]) > bigger:
				bigger = np.median(myhash2[pair])
			median_diff_perc = abs( np.median(myhash[pair])-np.median(myhash2[pair]) ) / bigger *100
			u, p_val = stat.mannwhitneyu(myhash[pair], myhash2[pair])
			outfile.write(str(j)+";"+str(k)+";"+str( p_val )+";"+str( median_diff )+";"+str( std_diff )+";"+str( median_diff_perc )+"\n")
outfile.close()

