#
#
import numpy as np;
import matplotlib.pyplot as plt;
import sys;

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


ens_type = "ensemble" # "ensemble" or "random_pool"
data_folder = "data/"+ens_type+"/"
res_folder = "results/"+ens_type+"/"
ensemble_file = "PED00001e002.pdb"
atomtype = "CA"
first_res = -1
last_res = 90

myhash={}
infile = open(data_folder+ensemble_file, 'r')
lines = []
for line0 in infile:
	if line0[0:4] == "ATOM" or line0[0:6] == "HETATM":
		if line0[13:16].strip() == 'CA':
			lines.append(line0)
infile.close()

print "Processed models:\n"
mod_num = 1
for line in lines:
	line = line.strip()
	j = int(line[22:26])
	if line[13:16].strip() == atomtype:
		lines2 = lines
		mod_num2 = 1
		for line2 in lines2:
			line2 = line2.strip()
			k = int(line2[22:26])
			if k > j:
				if line2[13:16].strip() == atomtype:
					if mod_num == mod_num2:
						mydist = calc_dist(line, line2)
						pair = (j,k)
						mylist = []
						if pair in myhash:
							mylist = myhash[pair]
						mylist.append(mydist)
						myhash[pair] = mylist
					if k == last_res:
						mod_num2 += 1
		if j == last_res:
			print mod_num
			mod_num += 1

outfile = open(res_folder+ensemble_file[0:-4]+"_"+atomtype+"_distance_distributions.txt", 'w')
matrix = np.empty((last_res-first_res+1, last_res-first_res+1))
matrix.fill(0.0)
for j in range(first_res, last_res):
	for k in range(j, last_res+1):
		if j == k:
			matrix[j-first_res][k-first_res] = 0.0
		else:
			pair = (j,k)
			if pair in myhash.keys():
				outfile.write(str(j)+";"+str(k)+";"+str(myhash[pair])+"\n")
				matrix[j-first_res][k-first_res] = np.median( myhash[pair] )
				matrix[k-first_res][j-first_res] = np.std( myhash[pair] )
			else:
				print j, k, "-> no distribution recorded"
outfile.close()
plt.imshow(matrix, cmap='CMRmap_r', interpolation='None', extent=[first_res,last_res,last_res,first_res])
plt.suptitle(ensemble_file[0:-4]+" "+atomtype+" distance matrix with st.dev", fontsize=18)
plt.title('top: median,  bottom: st.dev', fontsize=12)
plt.xlabel('residue number')
plt.ylabel('residue number')
plt.colorbar()
plt.savefig(res_folder+ensemble_file[0:-4]+"_"+atomtype+"_ensemble_dist_heatmap_w_std.png", dpi=600)
plt.clf()
