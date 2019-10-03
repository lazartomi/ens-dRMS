#
#
import scipy.stats as stat
import numpy as np
import matplotlib.pyplot as plt


for n in range(1,6):
	
	myhash_median_diff = {}
	myhash_median_diff_perc = {}
	myhash_median_diff_signif = {}
	myhash_median_diff_perc_signif = {}
	myhash_std_diff = {}
	
	sum_d_squared = 0.0
	ensemble = "6AAC-"+str(n) #7AAC
	first_res = 1
	last_res = 130 #130 for 6AAC, 132 for 7AAC

	infilename = "../PED6AAC-k18_pool_files/"+ensemble+"-pool_distance_distribution_differences.txt" #
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
			if p_val <= 0.05:
				sum_d_squared += np.square( median_diff )
	infile.close()
	N = last_res+1-first_res
	Npairs = N*(N-1)/2
	ensRMSD = np.sqrt(sum_d_squared / Npairs)
	print ensRMSD


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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=10, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool difference matrix with median & st.dev", fontsize=18) #PED7AAC-ntail_pool
	plt.title('top: d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED6AAC-k18_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_std_wo_test.png", dpi=600) ## PED7AAC-ntail_pool_files/PED7AAC-ntail_pool...
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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=20, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool normalized difference matrix", fontsize=18) #PED6AAC-k18_pool
	plt.title('top: %d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED6AAC-k18_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_perc_std_wo_test.png", dpi=600) # PED6AAC-k18_pool_files/PED6AAC-k18_pool...
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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=10, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool difference matrix with median & st.dev", fontsize=18) #PED7AAC-ntail_pool
	plt.title('top: d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED6AAC-k18_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_std.png", dpi=600) ## PED7AAC-ntail_pool_files/PED7AAC-ntail_pool...
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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=20, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool normalized difference matrix", fontsize=18) #PED6AAC-k18_pool
	plt.title('top: %d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED6AAC-k18_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_perc_std.png", dpi=600) # PED6AAC-k18_pool_files/PED6AAC-k18_pool...
	plt.clf()
	
	
#--------------------------


for n in range(1,6):
	
	myhash_median_diff = {}
	myhash_median_diff_perc = {}
	myhash_median_diff_signif = {}
	myhash_median_diff_perc_signif = {}
	myhash_std_diff = {}
	
	sum_d_squared = 0.0
	ensemble = "7AAC-"+str(n) #7AAC
	first_res = 1
	last_res = 132 #130 for 6AAC, 132 for 7AAC

	infilename = "../PED7AAC-ntail_pool_files/"+ensemble+"-pool_distance_distribution_differences.txt" #
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
			if p_val <= 0.05:
				sum_d_squared += np.square( median_diff )
	infile.close()
	N = last_res+1-first_res
	Npairs = N*(N-1)/2
	ensRMSD = np.sqrt(sum_d_squared / Npairs)
	print ensRMSD


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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=10, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool difference matrix with median & st.dev", fontsize=18) #PED7AAC-ntail_pool
	plt.title('top: d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED7AAC-ntail_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_std_wo_test.png", dpi=600) ## PED7AAC-ntail_pool_files/PED7AAC-ntail_pool...
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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=20, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool normalized difference matrix", fontsize=18) #PED6AAC-k18_pool
	plt.title('top: %d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED7AAC-ntail_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_perc_std_wo_test.png", dpi=600) # PED6AAC-k18_pool_files/PED6AAC-k18_pool...
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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=10, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool difference matrix with median & st.dev", fontsize=18) #PED7AAC-ntail_pool
	plt.title('top: d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED7AAC-ntail_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_std.png", dpi=600) ## PED7AAC-ntail_pool_files/PED7AAC-ntail_pool...
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
	
	
	plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', vmax=20, extent=[first_res,last_res,last_res,first_res])
	plt.suptitle(ensemble+"-pool normalized difference matrix", fontsize=18) #PED6AAC-k18_pool
	plt.title('top: %d(median),  bottom: d(st.dev)', fontsize=12)
	plt.xlabel('residue number')
	plt.ylabel('residue number')
	plt.colorbar()
	plt.savefig("../PED7AAC-ntail_pool_files/"+ensemble+"-pool_ensemble_diff_heatmap_w_median_perc_std.png", dpi=600) # PED6AAC-k18_pool_files/PED6AAC-k18_pool...
	plt.clf()
	
	