# ens-dRMS

To address the challenge of comparing ensembles of disordered proteins, we developed superposition-independent measures for evaluating the local and global similarity between two ensembles. The local similarity between specific regions of the polypeptide is evaluated from the differences between the distance distributions of individual residue pairs, and their statistical significance.  The global similarity is captured by the ens-dRMS measure, an RMSD-like quantity representing the root mean square difference between the medians of the inter-residue distance distributions of the two ensembles.


The distance-based metrics are complemented with several classical measures applied to individual ensembles. The local backbone variability is quantified by the distributions of the average backbone RMSD values of 5-residue segments along the polypeptide, computed over pairs of conformations in each ensemble. Global conformational parameters are also quantified from the distribution the radius of gyration (Rg) of conformations in individual ensembles.

<b>The analysis pipeline is made available as a Google Colab notebook (scroll to the bottom) and as separate scripts for local execution.</b>

<u>Desciption of the scripts (recommended to execute in the following order):</u>

1. dist_distrib_CA_fast.py

  Calculates the residue-resiude distance distributions for each amino acid pairs. Input variables include the data folder, the type (ensemble/pool) the output folder, and the first and last amino acids of the ensemble construct. By default the amino acid is represented by its CA atom. It outputs raw data in a file in the results folder, and also outputs the progress to the terminal.
  
2.  ens-RMSD.py

  Calculates the difference between residue-resiude distance distributions. It requires the distance distribution files output by the previous script. Other inputs besides the folder and filename, include the first and last amino acids of the ensemble construct. By default the amino acid is represented by its CA atom – obviously, it should match with the setup of the previous script. This script outputs files with raw data in the results folder (differences between distance distributions).
  
3. plot_diff_ens-RMSD.py

  In this scipt, inputs besides the folder and filename include the first and last amino acids of the ensemble construct. The script prints the ens-dRMS values (former versions termed it ensRMSD) to the terminal (best be redirected to a file) and plots the difference maps in the results folder. Multiple versions of the calculation will be carried out: with or without the normalization (percentages instead of absolute values), moreover with or without masking out the not significantly different amino acid pairs.

4. compare_Rg_distributions.py

  Compares the Rg distribution of two ensembles (one can be a random pool). Input variables include the type (ensemble/pool) the output folder, first and last residues. The script plots the corresponding Rg distributions.
  
5. local_superposition.py

  This script needs Biopython's Bio.PDB to be installed. Input variables are given in the script (data and results folder, type, first and last amino acid). It outputs the RMSD for pentapeptide regions of a sliding window in a data file created in the results folder.

6. plot_local_RMSDs.py

  It plots the RMSD for pentapeptide regions of a sliding window. For this, the script needs the output of the previous script as input (plus first and last amino acids of the construct). The output plot is a lineplot of the mean and the 95 percentiles with grey shade along the sequence.
  

# Google Colab Notebook:

ens_dRMS_Google_Colab.ipynb

  All the pipeline steps are realized in the notebook (with minimal changes) for accessibility and easy exploration of the pipeline elements.
  
