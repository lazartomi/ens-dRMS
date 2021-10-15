#
#
import scipy.stats as stat
import numpy as np
import matplotlib.pyplot as plt


ens_type = "ensemble" # "ensemble" or "random_pool"
res_folder = "results/"+ens_type+"/"
ensemble = "PED00001e001"
ensemble2 = "PED00001e002"
atomtype = "CA"
first_res = -1
last_res = 90


myhash_median_diff = {}
myhash_median_diff_perc = {}
myhash_median_diff_signif = {}
myhash_median_diff_perc_signif = {}
myhash_std_diff = {}

sum_d_squared = 0.0
sum_d_squared_all = 0.0
infilename = "results/"+ens_type+"/"+ensemble+"-"+ensemble2+"_"+atomtype+"_dist_distrib_differences.txt"
infile = open(infilename, 'r')
for line in infile:
	line = line.strip()
	if len(line)>1:
		array = line.split(";")
		pair = ( int(array[0]), int(array[1]) )
		p_val = float( array[2] )
		median_diff = float( array[3] )
		std_diff = float( array[4] )
		median_diff_perc = float( array[5] )
		#median
		mylist = []
		if pair in myhash_median_diff.keys():
			mylist = myhash_median_diff[pair]
		mylist.append(median_diff)
		myhash_median_diff[pair] = mylist
		#median_signif
		mylist = []
		if pair in myhash_median_diff_signif.keys():
			mylist = myhash_median_diff_signif[pair]
		if p_val <= 0.05:
			mylist.append(median_diff)
		else:
			mylist.append(0.0)
		myhash_median_diff_signif[pair] = mylist
		#std
		mylist = []
		if pair in myhash_std_diff.keys():
			mylist = myhash_std_diff[pair]
		mylist.append(std_diff)
		myhash_std_diff[pair] = mylist
		#median diff perc
		mylist = []
		if pair in myhash_median_diff_perc.keys():
			mylist = myhash_median_diff_perc[pair]
		mylist.append(median_diff_perc)
		myhash_median_diff_perc[pair] = mylist
		#median diff perc signif
		mylist = []
		if pair in myhash_median_diff_perc_signif.keys():
			mylist = myhash_median_diff_perc_signif[pair]
		if p_val <= 0.05:
			mylist.append(median_diff_perc)
		else:
			mylist.append(0.0)
		myhash_median_diff_perc_signif[pair] = mylist
		#ensRMSD
		sum_d_squared_all += np.square( median_diff )
		if p_val <= 0.05:
			sum_d_squared += np.square( median_diff )
infile.close()
N = last_res+1-first_res
Npairs = N*(N-1)/2
ensRMSD = np.sqrt(sum_d_squared / Npairs)
print "ensRMSD="+str(ensRMSD)+" A     (masked)"
ensRMSD_all = np.sqrt(sum_d_squared_all / Npairs)
print "ensRMSD="+str(ensRMSD_all)+" A     (unmasked - not recommended)"


matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			matrix[j-first_res][k-first_res] = np.mean(myhash_median_diff[pair] )
			matrix[k-first_res][j-first_res] = np.mean(myhash_std_diff[pair] ) 


plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', extent=[first_res,last_res,last_res,first_res])
plt.suptitle("Difference matrix with median & st.dev", fontsize=18)
plt.title('top: d(median),  bottom: d(st.dev)', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig("results/"+ens_type+"/"+ensemble+"-"+ensemble2+"_"+atomtype+"_diff_heatmap_w_median_std_wo_test.png", dpi=600) ## PED7AAC-ntail_pool_files/PED7AAC-ntail_pool...
plt.clf()


matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			matrix[j-first_res][k-first_res] = np.mean(myhash_median_diff_perc[pair] )
			matrix[k-first_res][j-first_res] = np.mean(myhash_std_diff[pair] ) 


plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', extent=[first_res,last_res,last_res,first_res])
plt.suptitle("Normalized difference matrix", fontsize=18) 
plt.title('top: %d(median),  bottom: d(st.dev)', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig("results/"+ens_type+"/"+ensemble+"-"+ensemble2+"_"+atomtype+"_diff_heatmap_w_median_perc_std_wo_test.png", dpi=600) # PED6AAC-k18_pool_files/PED6AAC-k18_pool...
plt.clf()

matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			matrix[j-first_res][k-first_res] = np.mean(myhash_median_diff_signif[pair] )
			matrix[k-first_res][j-first_res] = np.mean(myhash_std_diff[pair] ) 


plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', extent=[first_res,last_res,last_res,first_res])
plt.suptitle("Difference matrix with median & st.dev - Masked", fontsize=18)
plt.title('top: d(median),  bottom: d(st.dev)', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig("results/"+ens_type+"/"+ensemble+"-"+ensemble2+"_"+atomtype+"_diff_heatmap_w_median_std.png", dpi=600) ## PED7AAC-ntail_pool_files/PED7AAC-ntail_pool...
plt.clf()


matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			matrix[j-first_res][k-first_res] = np.mean(myhash_median_diff_perc_signif[pair] )
			matrix[k-first_res][j-first_res] = np.mean(myhash_std_diff[pair] ) 


plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', extent=[first_res,last_res,last_res,first_res])
plt.suptitle("Normalized difference matrix - Masked", fontsize=18)
plt.title('top: %d(median),  bottom: d(st.dev)', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig("results/"+ens_type+"/"+ensemble+"-"+ensemble2+"_"+atomtype+"_diff_heatmap_w_median_perc_std.png", dpi=600) # PED6AAC-k18_pool_files/PED6AAC-k18_pool...
plt.clf()
	