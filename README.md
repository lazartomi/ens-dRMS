# ens-dRMS


To address the challenge of comparing ensembles of disordered proteins, we developed superposition-independent measures for evaluating the local and global similarity between two ensembles. The local similarity between specific regions of the polypeptide is evaluated from the differences between the distance distributions of individual residue pairs, and their statistical significance.  The global similarity is captured by the ens-dRMS measure, an RMSD-like quantity representing the root mean square difference between the medians of the inter-residue distance distributions of the two ensembles.


The distance-based metrics are complemented with several classical measures applied to individual ensembles.  The local backbone variability is quantified by the distributions of the average backbone RMSD values of 5-residue segments along the polypeptide, computed over pairs of conformations in each ensemble. Local conformational preferences are evaluated by the frequencies of backbone (phi, psi torsion angles of individual residues, mapped onto regions of the Ramachandran map corresponding to secondary structure motifs. Global conformational parameters are also quantified from the distribution the radius of gyration (Rg) of conformations in individual ensembles.
