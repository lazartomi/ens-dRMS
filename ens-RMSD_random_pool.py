#
#
import scipy.stats as stat
import numpy as np
'''
for n in range(1,100):#100
	for m in range(n+1, 101):#101
		ensemble = "PED6AAC-k18_pool-"+str(n)
		ensemble2 = "PED6AAC-k18_pool-"+str(m)
		first_res = 1
		last_res = 130

		myhash = {}
		infile = open("PED6AAC-k18_pool_files/"+ensemble+"_distance_distributions.txt", 'r')
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

		myhash2 = {}
		infile2 = open("PED6AAC-k18_pool_files/"+ensemble2+"_distance_distributions.txt", 'r')
		for line2 in infile2:
			line2 = line2.strip()
			if len(line2)>1:
				array2 = line2.split(";")
				pair2 = ( int(array2[0]), int(array2[1]) )
				numbers2 = array2[2][1:-1]
				number_array2 = numbers2.split(", ")
				num_array2 = []
				for i in number_array2:
					num_array2.append(float(i))
				myhash2[pair2] = num_array2
		infile2.close()

		outfilename = ""
		if m<10:
			outfilename = "PED6AAC-k18_pool_files/Distribution_differences/"+ensemble+"-"+ensemble2[-1:]+"_distance_distribution_differences.txt"
		elif m<100:
			outfilename = "PED6AAC-k18_pool_files/Distribution_differences/"+ensemble+"-"+ensemble2[-2:]+"_distance_distribution_differences.txt"
		else:
			outfilename = "PED6AAC-k18_pool_files/Distribution_differences/"+ensemble+"-"+ensemble2[-3:]+"_distance_distribution_differences.txt"
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
					#print
					outfile.write(str(j)+";"+str(k)+";"+str( p_val )+";"+str( median_diff )+";"+str( std_diff )+";"+str( median_diff_perc )+"\n")
		outfile.close()



'''
for n in range(1,100):#100
	for m in range(n+1, 101):#101
		ensemble = "PED7AAC-ntail_pool-"+str(n)
		ensemble2 = "PED7AAC-ntail_pool-"+str(m)
		first_res = 1
		last_res = 132

		myhash = {}
		infile = open("PED7AAC-ntail_pool_files/"+ensemble+"_distance_distributions.txt", 'r')
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

		myhash2 = {}
		infile2 = open("PED7AAC-ntail_pool_files/"+ensemble2+"_distance_distributions.txt", 'r')
		for line2 in infile2:
			line2 = line2.strip()
			if len(line2)>1:
				array2 = line2.split(";")
				pair2 = ( int(array2[0]), int(array2[1]) )
				numbers2 = array2[2][1:-1]
				number_array2 = numbers2.split(", ")
				num_array2 = []
				for i in number_array2:
					num_array2.append(float(i))
				myhash2[pair2] = num_array2
		infile2.close()

		outfilename = ""
		if m<10:
			outfilename = "PED7AAC-ntail_pool_files/Distribution_differences/"+ensemble+"-"+ensemble2[-1:]+"_distance_distribution_differences.txt"
		elif m<100:
			outfilename = "PED7AAC-ntail_pool_files/Distribution_differences/"+ensemble+"-"+ensemble2[-2:]+"_distance_distribution_differences.txt"
		else:
			outfilename = "PED7AAC-ntail_pool_files/Distribution_differences/"+ensemble+"-"+ensemble2[-3:]+"_distance_distribution_differences.txt"
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
					#print
					outfile.write(str(j)+";"+str(k)+";"+str( p_val )+";"+str( median_diff )+";"+str( std_diff )+";"+str( median_diff_perc )+"\n")
		outfile.close()

