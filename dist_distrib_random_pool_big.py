#
#
import numpy as np;
import matplotlib.pyplot as plt;

def calc_dist(line, line2):
	d = 0.0
	x1 = float(line[29:38])
	x2 = float(line2[29:38])
	y1 = float(line[38:46])
	y2 = float(line2[38:46])
	z1 = float(line[46:54])
	z2 = float(line2[46:54])

	d = np.sqrt( np.square(x1-x2)+np.square(y1-y2)+np.square(z1-z2) )
	return d;



folder = "../PED7AAC-ntail_pool/"
first_res = 1
last_res = 132
chain = " "

myhash = {}
for i in range(0, 20000):
	if i%100 == 0:
		print folder, i
	for j in range(first_res, last_res):
		infile = open(folder+str(i)+"a_132.pdb", 'r') #a_130.pdb
		for line in infile:
			line = line.strip()
			if len(line)>1:
				if line[0:4] == "ATOM" or line[0:6] == "HETATM":
					if line[21:22] == chain and line[13:16].strip() == "CA":
						if int(line[22:26]) == j:
							for k in range(j+1, last_res+1):
								infile2 = open(folder+str(i)+"a_132.pdb", 'r') #a_130.pdb
								for line2 in infile2:
									line2 = line2.strip()
									if len(line2)>1:
										if line2[0:4] == "ATOM" or line2[0:6] == "HETATM":
											if line2[21:22] == chain and line2[13:16].strip() == "CA":
												if int(line2[22:26]) == k:
													mydist = calc_dist(line, line2)
													pair = (j,k)
													mylist = []
													if pair in myhash:
														mylist = myhash[pair]
													mylist.append(mydist)
													myhash[pair] = mylist
													break;
								infile2.close()
		infile.close()

outfile = open("../PED7AAC-ntail_pool_files/"+folder[3:19]+"_distance_distributions.txt", 'w')
matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			if pair in myhash.keys():
				print j, k, np.median( myhash[pair] ), np.std( myhash[pair] )
				outfile.write(str(j)+";"+str(k)+";"+str(myhash[pair])+"\n")
				matrix[j-first_res][k-first_res] = np.median( myhash[pair] )
				matrix[k-first_res][j-first_res] = np.std( myhash[pair] )
			else:
				print j, k, "-> no distribution recorded"
outfile.close()
plt.imshow(matrix, cmap='gray', interpolation='None')
plt.suptitle(folder[3:19]+"- distance matrix with st.dev", fontsize=18)
plt.title('top: median,  bottom: st.dev', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig("../PED7AAC-ntail_pool_files/"+folder[3:19]+"_ensemble_dist_heatmap_w_std.png", dpi=600)
plt.clf()


folder = "../PED6AAC-k18_pool/"
first_res = 1
last_res = 130
chain = " "

myhash = {}
for i in range(0, 20000):
	if i%100 == 0:
		print folder, i
	for j in range(first_res, last_res):
		infile = open(folder+str(i)+"a_130.pdb", 'r') #a_130.pdb
		for line in infile:
			line = line.strip()
			if len(line)>1:
				if line[0:4] == "ATOM" or line[0:6] == "HETATM":
					if line[21:22] == chain and line[13:16].strip() == "CA":
						if int(line[22:26]) == j:
							for k in range(j+1, last_res+1):
								infile2 = open(folder+str(i)+"a_130.pdb", 'r') #a_130.pdb
								for line2 in infile2:
									line2 = line2.strip()
									if len(line2)>1:
										if line2[0:4] == "ATOM" or line2[0:6] == "HETATM":
											if line2[21:22] == chain and line2[13:16].strip() == "CA":
												if int(line2[22:26]) == k:
													mydist = calc_dist(line, line2)
													pair = (j,k)
													mylist = []
													if pair in myhash:
														mylist = myhash[pair]
													mylist.append(mydist)
													myhash[pair] = mylist
													break;
								infile2.close()
		infile.close()

outfile = open("../PED6AAC-k18_pool_files/"+folder[3:19]+"_distance_distributions.txt", 'w')
matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			if pair in myhash.keys():
				print j, k, np.median( myhash[pair] ), np.std( myhash[pair] )
				outfile.write(str(j)+";"+str(k)+";"+str(myhash[pair])+"\n")
				matrix[j-first_res][k-first_res] = np.median( myhash[pair] )
				matrix[k-first_res][j-first_res] = np.std( myhash[pair] )
			else:
				print j, k, "-> no distribution recorded"
outfile.close()
plt.imshow(matrix, cmap='gray', interpolation='None')
plt.suptitle(folder[3:19]+"- distance matrix with st.dev", fontsize=18)
plt.title('top: median,  bottom: st.dev', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig("../PED6AAC-k18_pool_files/"+folder[3:19]+"_ensemble_dist_heatmap_w_std2.png", dpi=600)
plt.clf()