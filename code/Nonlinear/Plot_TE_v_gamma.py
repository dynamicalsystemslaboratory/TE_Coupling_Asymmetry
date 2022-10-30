#from matplotlib.lines import _LineStyle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size" : 17,
    "lines.linewidth" : 2.0,
    "font.family": "arial"
    })

fdir = './figures_bin/'; 
fdir2 = './files_bin/'; 

a=0.3; b=0.4; eps=0.1; alph=2.0
r = 0.1; rstr = '0_1'; frstr = '_r'+rstr

k=1; l=1
fstr = 'a0_3_b0_4_eps0_1_alph2_sig0_1_Nt1mil'

gamarr = np.array([0.0,0.01,0.02,0.05]) #a=0.3,b=0.4
gamstr = np.array(['0','0_01','0_02','0_05'])

Ngam = np.size(gamarr)

netTExy = np.zeros(Ngam)
netTECxy = np.zeros(Ngam)
netTEC2xy = np.zeros(Ngam)

for igam in range(Ngam):
    gam = gamarr[igam]
    fstr2= '_gam' + gamstr[igam]

    str3= fdir2 + 'TE_bin_all_vs_r_' + fstr + fstr2 + '.csv' # actual dist
    DF = pd.read_csv(str3, index_col=0)
    DF3 = DF[(DF['r']==r)] 
    TExy = np.array(DF3['TExy']) 
    TEyx = np.array(DF3['TEyx']) 
    TECxy = np.array(DF3['TECxy']) 
    TECyx = np.array(DF3['TECyx']) 
    TEC2xy = np.array(DF3['TEC2xy']) 
    TEC2yx = np.array(DF3['TEC2yx']) 

    netTExy[igam] = TExy - TEyx
    netTECxy[igam] = TECxy - TECyx
    netTEC2xy[igam] = TEC2xy - TEC2yx

## Plot cond entropy for using X as target
plt.figure()
fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=(5, 3.5))
ax1.plot(gamarr,netTExy, linewidth=2, marker='*', markersize=10, label=r'$\mathrm{net \; TE}_{X → Y}$')
ax1.plot(gamarr,netTECxy, linewidth=2, marker='*', markersize=10, label=r'$\mathrm{net \; TE}^\mathrm{C1}_{X → Y}$')
ax1.plot(gamarr,netTEC2xy, linewidth=2, marker='*', markersize=10, label=r'$\mathrm{net \; TE}^\mathrm{C2}_{X → Y}$')
ax1.plot(gamarr,np.zeros(gamarr.shape), color='black',linewidth=2, linestyle='dashed')
ax1.set_xlabel(r'$\gamma$',fontsize=15)
ax1.set_ylabel(r'$\mathrm{net \; TE}_{X → Y},\mathrm{net \; TE}^\mathrm{C1}_{X → Y},\mathrm{net \; TE}^\mathrm{C2}_{X → Y}$',fontsize=13)
ax1.set_ylim([-0.01,0.03]) 
plt.legend(fontsize=13,loc='best', bbox_to_anchor=(0.5, 0.37, 0.5, 0.5))  
fig.savefig(fdir + 'netTE_all_vs_gamma_' + fstr + frstr + '.png', dpi=150, bbox_inches="tight")
fig.savefig(fdir + 'netTE_all_vs_gamma_' + fstr + frstr + '.eps', dpi=150, bbox_inches="tight")

