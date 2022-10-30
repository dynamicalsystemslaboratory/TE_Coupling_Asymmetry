#from matplotlib.lines import _LineStyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import time

from pylab import plot, show, savefig, xlim, figure, \
                  ylim, legend, boxplot, setp, axes
                  
plt.rcParams.update({
    "font.size" : 17,
    "lines.linewidth" : 2.0,
    "font.family": "arial"
    })

# set parameters
simN = 1
N = 20 # no. of simulations 

r=0.5;  rstr = '0_5' # for filenames
eps=0.1; k=1; l=1

fdir = './figures_bin/'; 
fdir2 = './files_bin_sims2/'; 
fstr = 'TE_vs_alph_eps0_1_sig0_1_r' + rstr + '_Nt10mil'

alpharr = np.array([0,0.5,1,2,4,6,8])
Nalph = np.size(alpharr)

fsims=glob.glob(fdir2+fstr+'*')
Nsims=len(fsims)
print(Nsims)

netTExy = np.zeros((Nalph,Nsims))
netTECxy = np.zeros((Nalph,Nsims))
netTEC2xy = np.zeros((Nalph,Nsims))

for ksim in range(0,Nsims):
    fname = fsims[ksim]

    DF = pd.read_csv(fname, index_col=0)
    alpharr = np.array(DF['alpha']) 
    TExy = np.array(DF['TExy']) 
    TEyx = np.array(DF['TEyx']) 
    TECxy = np.array(DF['TECxy']) 
    TECyx = np.array(DF['TECyx'])
    TEC2xy = np.array(DF['TEC2xy']) 
    TEC2yx = np.array(DF['TEC2yx'])

    netTExy[:,ksim] = TExy-TEyx
    netTECxy[:,ksim] = TECxy-TECyx
    netTEC2xy[:,ksim] = TEC2xy-TEC2yx

countTE = np.zeros(Nalph); countTEC = np.zeros(Nalph); countTEC2 = np.zeros(Nalph)

for ialph in range(0,Nalph):
    alph = alpharr[ialph]
    print(" alpha = ", alph, "---------------------")

    countTE[ialph] = np.sum((netTExy[ialph,:]>0)*1)
    countTEC[ialph] = np.sum((netTECxy[ialph,:]>0)*1)
    countTEC2[ialph] = np.sum((netTEC2xy[ialph,:]>0)*1)


## Plot cond entropy for using X as target
fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=(5, 3.5))
ax1.plot(alpharr,countTE*100.0/Nsims, linewidth=2, marker='*', markersize=10, label=r'$\mathrm{net \; TE}_{X → Y}$')
ax1.plot(alpharr,countTEC*100.0/Nsims, linewidth=2, marker='*', markersize=10, label=r'$\mathrm{net \; TE}^\mathrm{C1}_{X → Y}$')
ax1.plot(alpharr,countTEC2*100.0/Nsims, linewidth=2, marker='*', markersize=10, label=r'$\mathrm{net \; TE}^\mathrm{C2}_{X → Y}$')
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$Correct \; coupling \; (\%)$',fontsize=15)
ax1.set_ylim([0,105])
plt.legend(fontsize=12)
fig.savefig(fdir + 'percC_netTE_all3_vs_alph_' + fstr + '.eps', dpi=150, bbox_inches="tight")
fig.savefig(fdir + 'percC_netTE_all3_vs_alph_' + fstr + '.png')

