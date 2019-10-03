#
#
import Bio.PDB
import matplotlib.pyplot as plt
import sys
import pickle
import os
import numpy as np

infilename = "../../DISICL_analysis/PED6AAC-pdb/6AAC_DISICL.pdb1"

structure = Bio.PDB.PDBParser(QUIET=True).get_structure("E1-5", infilename)
print "Ensemble parsed"

aafrom = 2
aato = 130
phis = []
psis = []
aa_phis = {}
aa_psis = {}

for model in structure:
	for chain in model :
		poly = Bio.PDB.Polypeptide.Polypeptide(chain)
		num = aafrom
		for pair in poly.get_phi_psi_list():
			if None not in pair:
				phis.append(pair[0]*57.2958)
				psis.append(pair[1]*57.2958)
				if num not in aa_phis.keys():
					aa_phis[num] = []
					aa_psis[num] = []
				aa_phis[num].append(pair[0]*57.2958)
				aa_psis[num].append(pair[1]*57.2958)
				num += 1
'''
fig = plt.figure(frameon = False)
fig.set_size_inches(6,6)
plt.title(i)
plt.scatter(phis, psis, c='k', marker='.', alpha=0.2)
plt.xlim(-185,185)
plt.ylim(-185,185)
plt.xlabel("Phi")
plt.ylabel("Psi")
plt.savefig(i+"_rama.png", dpi=600)
plt.clf()
'''
frac_BetaStrand = []
frac_BetaPro = [] #PolyProII
frac_Alpha = []
frac_LeftAlpha = []
for k in range(aafrom,aato):
	betaStrand = 0.0
	betaPro = 0.0
	alpha = 0.0
	leftAlpha = 0.0
	for j in range(0,len(aa_phis[k]) ):
		if aa_phis[k][j] >= 0.0:
			leftAlpha += 1.0
		else:
			if aa_psis[k][j] < 50 and aa_psis[k][j] > -120:
				alpha += 1.0
			else:
				if aa_phis[k][j] < -100:
					betaStrand += 1.0
				else:
					betaPro += 1.0
	frac_BetaStrand.append( betaStrand / len(aa_phis[k]) )
	frac_BetaPro.append( betaPro / len(aa_phis[k]) )
	frac_Alpha.append( alpha / len(aa_phis[k]) )
	frac_LeftAlpha.append( leftAlpha / len(aa_phis[k]) )
#frac_BetaStrand.append( 0.25 )
#frac_BetaPro.append( 0.25 )
#frac_Alpha.append( 0.25 )
#frac_LeftAlpha.append( 0.25 )


'''

#for i in range(0,99):
	aafrom = 23
	aato = 45
	if i == 'Absinth':
		aafrom = 2
		aato =  24
	r_aa_phis = {}
	r_aa_psis = {}
	fig = plt.figure(frameon = False)
	fig.set_size_inches(6,6)
	for n in range(0+(200*i), 200+(200*i) ):
		phis = []
		psis = []
		infilename = "PED"+pool+"/"+str(n)+"a_"+str(aa_length)+".pdb" #6AAC: a_130.pdb, 7AAC: a_132.pdb
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure("PED"+pool+"-"+str(n), infilename)
		for model in structure:
			for chain in model :
				poly = Bio.PDB.Polypeptide.Polypeptide(chain)
				num = 2
				for pair in poly.get_phi_psi_list():
					if None not in pair:
						phis.append(pair[0]*57.2958)
						psis.append(pair[1]*57.2958)
						if num not in r_aa_phis.keys():
							r_aa_phis[num] = []
							r_aa_psis[num] = []
						r_aa_phis[num].append(pair[0]*57.2958)
						r_aa_psis[num].append(pair[1]*57.2958)
						num += 1
		
		plt.rcParams.update({'figure.max_open_warning': 0})
 		plt.title("PED"+pool+"-"+str(i+1))
		plt.scatter(phis, psis, c='k', marker='.', alpha=0.2)
		plt.xlim(-185,185)
		plt.ylim(-185,185)
		plt.xlabel("Phi")
		plt.ylabel("Psi")
		if n == (199+(200*i)):
			plt.savefig("PED"+pool+"_files/PED"+pool+"-"+str(i+1)+"_rama.png", dpi=600)
	plt.clf()
	r_frac_BetaStrand = [0.25]
	r_frac_BetaPro = [0.25] #PolyProII
	r_frac_Alpha = [0.25]
	r_frac_LeftAlpha = [0.25]
	for k in range(aafrom,aato):
		betaStrand = 0.0
		betaPro = 0.0
		alpha = 0.0
		leftAlpha = 0.0
		for j in range(0,len(r_aa_phis[k]) ):
			if r_aa_phis[k][j] >= 0.0:
				leftAlpha += 1
			else:
				if r_aa_psis[k][j] < 50 and r_aa_psis[k][j] > -120:
					alpha += 1
				else:
					if r_aa_phis[k][j] < -100:
						betaStrand += 1
					else:
						betaPro += 1
		r_frac_BetaStrand.append( betaStrand / len(r_aa_phis[k]) )
		r_frac_BetaPro.append( betaPro / len(r_aa_phis[k]) )
		r_frac_Alpha.append( alpha / len(r_aa_phis[k]) )
		r_frac_LeftAlpha.append( leftAlpha / len(r_aa_phis[k]) )
	r_frac_BetaStrand.append( 0.25 )
	r_frac_BetaPro.append( 0.25 )
	r_frac_Alpha.append( 0.25 )
	r_frac_LeftAlpha.append( 0.25 )
'''
	
'''
	fig = plt.figure(frameon = False)
	fig.set_size_inches(8,4.5)
	plt.title(i+" Beta")
	#plt.plot(r_frac_BetaStrand, 'k-')
	plt.plot(frac_BetaStrand, 'r-')	
	plt.xlim(0,23)
	plt.ylim(0.0,1.0)
	plt.xlabel("residues")
	plt.ylabel("fraction")
	plt.savefig(i+"_sec_struct_BS.png", dpi=600)
	plt.clf()

	plt.title(i+" PPII")
	#plt.plot(r_frac_BetaPro, 'k-')
	plt.plot(frac_BetaPro, 'r-')
	plt.xlim(0,23)
	plt.ylim(0.0,1.0)
	plt.xlabel("residues")
	plt.ylabel("fraction")
	plt.savefig(i+"_sec_struct_BP.png", dpi=600)
	plt.clf()

	plt.title(i+" Alpha")
	#plt.plot(r_frac_Alpha, 'k-')
	plt.plot(frac_Alpha, 'r-')
	plt.xlim(0,23)
	plt.ylim(0.0,1.0)
	plt.xlabel("residues")
	plt.ylabel("fraction")
	plt.savefig(i+"_sec_struct_aR.png", dpi=600)
	plt.clf()

	plt.title(i+" L-Alpha")
	#plt.plot(r_frac_LeftAlpha, 'k-')
	plt.plot(frac_LeftAlpha, 'r-')
	plt.xlim(0,23)
	plt.ylim(0.0,1.0)
	plt.xlabel("residues")
	plt.ylabel("fraction")
	plt.savefig(i+"_sec_struct_aL.png", dpi=600)
	plt.clf()
'''


ind = np.arange(aafrom,aato,1)
p1 = plt.bar(ind, frac_Alpha, 1.0, color='red')
p2 = plt.bar(ind, frac_LeftAlpha, 1.0, color = 'yellow', bottom=frac_Alpha)
p3 = plt.bar(ind, frac_BetaStrand, 1.0, color='blue', bottom=[a+b for a,b in zip(frac_Alpha, frac_LeftAlpha)] )
p4 = plt.bar(ind, frac_BetaPro, 1.0, color='green', bottom=[a+b+c for a,b,c in zip(frac_Alpha, frac_LeftAlpha, frac_BetaStrand)])
plt.ylabel('Fraction of dihedral angle assignments')
plt.xlabel('residues')
plt.xlim(aafrom, aato)
plt.ylim(0,1)
plt.legend( (p1[0], p2[0], p3[0], p4[0]), ('Alpha', 'L-alpha', 'Beta-str.', 'PPII'),
		bbox_to_anchor=(0.115, 0.99, 0.8, 0.002), loc=8, ncol=4 )
plt.savefig("PED6AAC-12345_sec_struct_every.png", dpi=600)
plt.clf()

print "Figure plotted"

