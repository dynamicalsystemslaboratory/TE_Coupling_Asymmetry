# TE_Coupling_Asymmetry

This readme file will help replicate the results in the paper titled: "Constructive critique to use of transfer entropy in inference of bidirectional interactions in heterogeneous systems".

There are 2 main directories:

Data: Contains links to the sources of the two real world example datasets used in this work

Code: Contains all the codes required to replicate the results of this work.


Code directory contains the following sub-directories, each for a different problem (analytical, numerical, or data-based):

1. Bivariate_Linear:

- “Linear_system_2D.m” is a code that calculates transfer entropy, controlled transfer entropy, and fully controlled transfer entropy for bivariate linear system at a range of epsilon & alpha values and compares with true dominant coupling direction to determine percentage of incorrect inferences for each measure
- “Linear_system_ND_Network.m” is a code that calculates the same as above but for a multi-node network system 
- “Plot_Percentage.m” is a code that plots the percentage of incorrect inferences vs. alpha
- “Plot_Percentage_Network.m” is a code that plots the percentage of incorrect inferences vs. alpha for multi-node network system
- “Plot_lineplots_Linear_system.m” is a code that calculates and plots transfer entropy, controlled transfer entropy, fully controlled transfer entropy, and self-regulation entropy for bivariate linear system at a given a,b,epsilon value for a range of alpha values
- “Plot_ab_contour_Linear_system_2D.m” is a code that calculates and plots contour plots of net transfer entropy, controlled transfer entropy, and fully controlled transfer entropy for bivariate linear system at a given epsilon & alpha value in the realizable (a,b) domain
- “Linear_system_Derivations.nb” is a Mathematica notebook that shows relevant derivations of bivariate/network linear system

2. Physiology:

- “TE_vs_r.py” is a code that computes transfer entropy, controlled transfer entropy, and fully controlled transfer entropy at different values of bandwidth r
- “plot_allTE.py” is a code that plots all measures as a function of bandwidth r
- “transfer_entropy_kde_all.py” is a function that computes transfer entropy, controlled transfer entropy, and fully controlled transfer entropy using kernel density estimation
- “compute_entropy.py” is a function that computes joint entropy of any given dimension

3. Pigeons:

- “WriteTE_all.py” is a code that computes conditional entropies, transfer entropy, controlled transfer entropy, and fully controlled transfer entropy for different generations and pairs of pigeons
- “cond_entropy_sym_k.py” is a function that computes conditional entropies
- “transfer_entropy_all_Hcond_sym_k.py” is a function that computes additional conditional entropies, transfer entropy, controlled transfer entropy, and fully controlled transfer entropy
- “TS_slope.py” is a code that creates the box plots of net transfer entropy, controlled transfer entropy, and fully controlled transfer entropy as a function of generations and computes the Theil-Sen slope and tests its significance
- “Variability.py” is a code that performs Mixed design ANOVA of mean deviation for different measures and also performs posthoc tests for pairwise comparison 
- “compute_entropy.py” is a function that computes joint entropy of any given dimension

4. Chaotic_Tent_Map:
 
- “TE2_v_alph.py” is a code that computes transfer entropy and general controlled transfer entropy for different k,l at different values of alpha
- “TE3_v_alph.py” is a code that computes transfer entropy, controlled transfer entropy, and fully controlled transfer entropy at different values of alpha
- “transfer_entropy2.py” is a function called to compute the value of  transfer entropy and general controlled transfer entropy
- “transfer_entropy3.py” is a function called to compute the value of  transfer entropy, controlled transfer entropy, and fully controlled transfer entropy
- “compute_entropy.py” is a function that computes joint entropy of any given dimension
- “Plot_netTE_v_alph.py” is a code that plots the net transfer entropy, controlled transfer entropy, and fully controlled transfer entropy at different values of alpha
- “Plot_perc_TE_v_alph.py” is a code that plots the percentage correct inference by net transfer entropy, controlled transfer entropy, and fully controlled transfer entropy at different values of alpha

5. Nonlinear:

- “TEall_vs_r_gam.py” is a code that computes transfer entropy, controlled transfer entropy, and fully controlled transfer entropy for different values of binwidth r and nonlinearity coefficient gamma
- “Plot_TE_v_gamma.py” is a code that plots transfer entropy, controlled transfer entropy, and fully controlled transfer entropy as a function of nonlinearity coefficient gamma
- “transfer_entropy3.py” is a function called to compute the value of  transfer entropy, controlled transfer entropy, and fully controlled transfer entropy
- “compute_entropy.py” is a function that computes joint entropy of any given dimension


NOTE: Please create the relevant files/ and figures/ folders within each of these subdirectories if executing the codes as is.
