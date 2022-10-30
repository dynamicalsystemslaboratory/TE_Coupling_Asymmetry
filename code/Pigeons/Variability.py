import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sc
from scipy import stats

from pylab import plot, show, savefig, xlim, figure, \
                  ylim, legend, boxplot, setp, axes

import pingouin as pg

lw=2.0
fsz = (9.5,4.5)

plt.rcParams.update({
    "font.size" : 22,
    "lines.linewidth" : lw,
    "font.family": "arial"
    })


dt = 0.2 # From paper 
k = 1
l = 1

minGen = 2; maxGen = 5
minRelease = 1; maxRelease = 12  # for each pair

# Modify string names for writing
fdir= './figures2/'; fstr = '_SS1_k'+str(k)+'_l'+str(l)+'.eps' 

# Sampling 
isamp = 0
samptarr = np.array([0.2]) # make sure sampt/dt = integer 
sampN = int(samptarr[isamp]/dt) # skip these many,i.e. downsampling fac

## Read actual MI,TE, ...
str3= './files/TE_all_Hcond_SS1_k'+str(k)+'_l'+str(l)+'.csv' # actual dist
DF3 = pd.read_csv(str3, index_col=0)
 
n1 = np.nan*np.ones(maxGen-minGen+1) 
SD_netTE_EN = np.nan*np.ones(maxGen-minGen+1)
SD_netTE_star_EN = np.nan*np.ones(maxGen-minGen+1)
SD_netTE_star2_EN = np.nan*np.ones(maxGen-minGen+1)
MAD_netTE_EN = np.nan*np.ones(maxGen-minGen+1)
MAD_netTE_star_EN = np.nan*np.ones(maxGen-minGen+1)
MAD_netTE_star2_EN = np.nan*np.ones(maxGen-minGen+1)
SDAD_netTE_EN = np.nan*np.ones(maxGen-minGen+1)
SDAD_netTE_star_EN = np.nan*np.ones(maxGen-minGen+1)
SDAD_netTE_star2_EN = np.nan*np.ones(maxGen-minGen+1)

Genarr = np.array([]) 
Nsubj = np.nan*np.ones(maxGen-minGen+1)
score1 = np.array([]); score2 = np.array([]); score3 = np.array([])

# Gen wise
for iGen in range(minGen,maxGen+1):

    print('Gen = ', iGen)

    # Get TE all for iGen
    df3 = DF3[(DF3['Gen']==iGen)]

    HcondE_gen = np.array(df3['HcondE']); HcondN_gen = np.array(df3['HcondN']); 

    TE_EN_gen = np.array(df3['TE_EN']); TE_NE_gen = np.array(df3['TE_NE'])
    TE_star_EN_gen = np.array(df3['TE_star_EN']); TE_star_NE_gen = np.array(df3['TE_star_NE'])
    TE_star2_EN_gen = np.array(df3['TE_star2_EN']); TE_star2_NE_gen = np.array(df3['TE_star2_NE'])
    HcondNE_gen = np.array(df3['HcondNE']); HcondEN_gen = np.array(df3['HcondEN']); 
    HcondNE2_gen = np.array(df3['HcondNE2']); HcondEN2_gen = np.array(df3['HcondEN2']); 

    n1[iGen-minGen] = TE_EN_gen.shape[0]; print('For this gen, n = ',n1[iGen-minGen])

    # normalized to compare variability
    netTE_EN_gen=(TE_EN_gen/HcondN_gen) - (TE_NE_gen/HcondE_gen)
    netTE_star_EN_gen=(TE_star_EN_gen/HcondNE_gen) - (TE_star_NE_gen/HcondEN_gen)
    netTE_star2_EN_gen=(TE_star2_EN_gen/HcondNE2_gen) - (TE_star2_NE_gen/HcondEN2_gen)

    print('norm TE = ',np.min((TE_EN_gen/HcondN_gen)),'-',np.max((TE_EN_gen/HcondN_gen)))
    print('norm TE* = ',np.min((TE_star_EN_gen/HcondNE_gen)),'-',np.max((TE_star_EN_gen/HcondNE_gen)))
    print('norm TE** = ',np.min((TE_star2_EN_gen/HcondNE2_gen)),'-',np.max((TE_star2_EN_gen/HcondNE2_gen)))

    # Std dev
    SD_netTE_EN[iGen-minGen] = np.std(netTE_EN_gen)
    SD_netTE_star_EN[iGen-minGen] = np.std(netTE_star_EN_gen)
    SD_netTE_star2_EN[iGen-minGen] = np.std(netTE_star2_EN_gen)

    # Mean absolute deviation
    MAD_netTE_EN[iGen-minGen] = np.mean( abs(netTE_EN_gen - np.mean(netTE_EN_gen)) )
    MAD_netTE_star_EN[iGen-minGen] = np.mean( abs(netTE_star_EN_gen - np.mean(netTE_star_EN_gen)) )
    MAD_netTE_star2_EN[iGen-minGen] = np.mean( abs(netTE_star2_EN_gen - np.mean(netTE_star2_EN_gen)) )

    # SD of absolute deviation
    SDAD_netTE_EN[iGen-minGen] = np.std( abs(netTE_EN_gen - np.mean(netTE_EN_gen)) )
    SDAD_netTE_star_EN[iGen-minGen] = np.std( abs(netTE_star_EN_gen - np.mean(netTE_star_EN_gen)) )
    SDAD_netTE_star2_EN[iGen-minGen] = np.std( abs(netTE_star2_EN_gen - np.mean(netTE_star2_EN_gen)) )

    # For mixed design ANOVA of abs difference from mean
    score1 = np.concatenate((score1,abs(netTE_EN_gen - np.mean(netTE_EN_gen))) ,axis=0)
    score2 = np.concatenate((score2,abs(netTE_star_EN_gen - np.mean(netTE_star_EN_gen))) ,axis=0)
    score3 = np.concatenate((score3,abs(netTE_star2_EN_gen - np.mean(netTE_star2_EN_gen))) ,axis=0)
    Nsubj[iGen-minGen] = netTE_EN_gen.shape[0]

    Genarr = np.concatenate((Genarr,iGen*np.ones(TE_EN_gen.shape)),axis=0)
    print('Gen ',iGen,' -----> ')

    stat, p = sc.stats.levene(netTE_EN_gen, netTE_star_EN_gen, netTE_star2_EN_gen, center='mean')
    print('Levenes test of x : p-value = ',p , ', F stat = ', stat)

    stat, p = sc.stats.f_oneway(abs(netTE_EN_gen - np.mean(netTE_EN_gen)), \
                    abs(netTE_star_EN_gen - np.mean(netTE_star_EN_gen)), \
                    abs(netTE_star2_EN_gen - np.mean(netTE_star2_EN_gen)) )
    print('one-way ANOVA of abs dev of x : p-value = ',p , ', F stat = ', stat)

## Mixed design ANOVA of abs difference scores (Levene's test)
Ntime = maxGen-minGen+1
print('Score1 size = ', score1.shape, ', 2 = ', score2.shape, ', 3 = ', score3.shape)
print('Nsubj = ', Nsubj)

df = pd.DataFrame({'Scores': np.r_[score1, score2, score3],
                    'Time': np.r_[Genarr,Genarr,Genarr],
                    'Group': np.repeat(['TE','TE*','TE**'], np.sum(Nsubj)),
                    'Subject': np.r_[np.arange(Nsubj[0]),
                                     np.arange(Nsubj[1]),
                                     np.arange(Nsubj[2]),
                                     np.arange(Nsubj[3]),
                                     np.arange(Nsubj[0],Nsubj[0]+Nsubj[0]),
                                     np.arange(Nsubj[1],Nsubj[1]+Nsubj[1]),
                                     np.arange(Nsubj[2],Nsubj[2]+Nsubj[2]),
                                     np.arange(Nsubj[3],Nsubj[3]+Nsubj[3]),
                                     np.arange(Nsubj[0]+Nsubj[0],Nsubj[0]+Nsubj[0]+Nsubj[0]),
                                     np.arange(Nsubj[1]+Nsubj[1],Nsubj[1]+Nsubj[1]+Nsubj[1]),
                                     np.arange(Nsubj[2]+Nsubj[2],Nsubj[2]+Nsubj[2]+Nsubj[2]),
                                     np.arange(Nsubj[3]+Nsubj[3],Nsubj[3]+Nsubj[3]+Nsubj[3])] })

print(df.head())

# Compute the two-way mixed-design ANOVA
aov = pg.mixed_anova(dv='Scores', within='Time', between='Group', subject='Subject', data=df)
# Pretty printing of ANOVA summary
pg.print_table(aov)

# Posthocs test to get significance in different Generations
posthocs = pg.pairwise_ttests(dv='Scores', within='Time', between='Group',
                              subject='Subject', data=df)
pg.print_table(posthocs)

# bar plots in grouped manner of bar type
x = np.arange(4)
y1 = MAD_netTE_EN 
y2 = MAD_netTE_star_EN 
y3 = MAD_netTE_star2_EN 
y1std = SDAD_netTE_EN/np.sqrt(n1) 
y2std = SDAD_netTE_star_EN/np.sqrt(n1) 
y3std = SDAD_netTE_star2_EN/np.sqrt(n1) 

width = 0.2
plt.figure()
ticks = ['2', '3', '4', '5']
plt.bar(x-0.2, y1, width, color='red', yerr = y1std, align='center', alpha=0.7, ecolor='darkred', capsize=4)
plt.bar(x, y2, width, color='blue', yerr = y2std, align='center', alpha=0.7, ecolor='mediumblue', capsize=4)
plt.bar(x+0.2, y3, width, color='green', yerr = y3std, align='center', alpha=0.7, ecolor='darkgreen', capsize=4)
plt.ylim([0,0.018]); plt.yticks(np.arange(0,0.021,0.003))
plt.xticks(np.arange(0, len(ticks)), ticks)
plt.xlabel('Generation'); plt.ylabel('mean absolute deviation')
plt.legend([r"$\mathrm{net \; TE}^{'}_{E → N}$",\
            r"$\mathrm{net \; TE}^\mathrm{C1'}_{E → N}$",\
            r"$\mathrm{net \; TE}^\mathrm{C2'}_{E → N}$"],loc=(1.02,0.35))
plt.savefig(fdir + 'barplot_MAD_mean_netTE_TEC_TEFC'+ fstr, bbox_inches='tight')
# plt.show()
