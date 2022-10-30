#from matplotlib.lines import _LineStyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pylab import plot, show, savefig, xlim, figure, \
                  ylim, legend, boxplot, setp, axes
                  

plt.rcParams.update({
    "font.size" : 17,
    "lines.linewidth" : 2.0,
    "font.family": "arial"
    })


fdir = './figures_bin/'; 
fdir2 = './files_bin2/'; 

r=0.5; str1 = '0_5' # for filenames

a=0.4; b=0.5; eps=0.1

k=1; l=1
fstr = 'a0_3_b0_4_eps0_1_sig0_1_r' + str1 + '_k' + str(k) + '_l' + str(l) + '_Nt10mil'

str3= fdir2 + 'TE_vs_alph_' + fstr + '.csv' # actual dist
DF3 = pd.read_csv(str3, index_col=0)
alpharr = np.array(DF3['alpha']) 
TExy = np.array(DF3['TExy']) 
TEyx = np.array(DF3['TEyx']) 
TECxy = np.array(DF3['TECxy']) 
TECyx = np.array(DF3['TECyx']) 

netTExy = TExy - TEyx
netTECxy = TECxy - TECyx

k=2; l=1
fstrB = 'a0_3_b0_4_eps0_1_sig0_1_r' + str1 + '_k' + str(k) + '_l' + str(l) + '_Nt10mil'

str3= fdir2 + 'TE_vs_alph_' + fstrB + '.csv' # actual dist
DF3 = pd.read_csv(str3, index_col=0)
# df = DF3[(DF3['l']==il+1)] 
alpharr = np.array(DF3['alpha']) 
TExy = np.array(DF3['TExy']) 
TEyx = np.array(DF3['TEyx']) 
TECxy = np.array(DF3['TECxy']) 
TECyx = np.array(DF3['TECyx']) 
Nalph = np.size(alpharr)

# netTExy = TExy - TEyx
netTEC2xy = TECxy - TECyx

## Plot 
fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=(5, 3.5))
ax1.plot(alpharr[0:Nalph-1],netTExy[0:Nalph-1], linewidth=2, label=r'$\mathrm{net \; TE}_{X → Y}$')  # r"$\mathrm{net \; TE^{C1'}}$"
ax1.plot(alpharr[0:Nalph-1],netTECxy[0:Nalph-1], linewidth=2, label=r'$\mathrm{net \; TE}^\mathrm{C1}_{X → Y}$')
ax1.plot(alpharr[0:Nalph-1],netTEC2xy[0:Nalph-1], linewidth=2, label=r'$\mathrm{net \; TE}^\mathrm{C2}_{X → Y}$')
ax1.plot(alpharr,np.zeros(alpharr.shape), color='black',linewidth=2, linestyle='dashed')
ax1.set_xlabel(r'$\alpha$')
ax1.set_ylabel(r'$\mathrm{net \; TE}_{X → Y},\mathrm{net \; TE}^\mathrm{C1}_{X → Y},\mathrm{net \; TE}^\mathrm{C2}_{X → Y}$',fontsize=13)
ax1.set_ylim([-0.3,0.1])
ax1.set_yticks([-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1])
plt.legend(fontsize=13)
fig.savefig(fdir + 'netTE_all3_vs_alph_' + fstr + '.eps', dpi=150, bbox_inches="tight")
fig.savefig(fdir + 'netTE_all3_vs_alph_' + fstr + '.png', dpi=150, bbox_inches="tight")

