import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sc
from scipy import stats

from pylab import plot, show, savefig, xlim, figure, \
                  ylim, legend, boxplot, setp, axes

lw=2.0
fsz = (7.2,6)

plt.rcParams.update({
    "font.size" : 24,
    "lines.linewidth" : lw,
    "font.family": "arial"
    })


def TS_slope(x,y):
    model = stats.mstats.theilslopes(y, x)
    print('Slope = ', model[0], 'Intercept = ', model[1])
    tau, p_value = stats.kendalltau(y, x)
    print('p_values: ')
    print('Kendal Tau test = ',p_value)
    test1 = stats.ranksums(x, y)
    print('Wilcoxon rank sums test = ', test1[1])
    test2 = stats.wilcoxon(x, y)
    print('Wilcoxon signed rank test = ', test2[1])

    # return in order - 
    R = np.array([model[0],model[1],p_value,test1[1],test2[1]])

    return R


dt = 0.2 # From paper 

minGen = 2; maxGen = 5
minRelease = 1; maxRelease = 12  # for each pair

# Modify string names for writing
fdir= './figures2/'; fstr = '_SS1_k1_l1.eps' 

# Sampling 
isamp = 0
samptarr = np.array([0.2]) # make sure sampt/dt = integer 
sampN = int(samptarr[isamp]/dt) # skip these many,i.e. downsampling fac

## Read actual MI,TE, ...
str3= './files/TE_all_Hcond_SS1_k1_l1.csv' # actual dist
DF3 = pd.read_csv(str3, index_col=0)
 
HcondE = np.array([]); HcondN = np.array([])
TE_EN = np.array([]); TE_NE = np.array([]) 
TE_star_EN = np.array([]); TE_star_NE = np.array([]) 
TE_star2_EN = np.array([]); TE_star2_NE = np.array([]) 
Genarr = np.array([]) 
HcondEN = np.array([]); HcondNE = np.array([])
HcondEN2 = np.array([]); HcondNE2 = np.array([])

HcondE_2D = np.nan*np.ones((maxGen-minGen+1,10*12)); HcondN_2D = np.nan*np.ones((maxGen-minGen+1,10*12)) 
TE_EN_2D = np.nan*np.ones((maxGen-minGen+1,10*12)); TE_NE_2D = np.nan*np.ones((maxGen-minGen+1,10*12)) 
TE_star_EN_2D = np.nan*np.ones((maxGen-minGen+1,10*12)); TE_star_NE_2D = np.nan*np.ones((maxGen-minGen+1,10*12)) 
TE_star2_EN_2D = np.nan*np.ones((maxGen-minGen+1,10*12)); TE_star2_NE_2D = np.nan*np.ones((maxGen-minGen+1,10*12)) 
HcondEN_2D = np.nan*np.ones((maxGen-minGen+1,10*12)); HcondNE_2D = np.nan*np.ones((maxGen-minGen+1,10*12)) 
HcondEN2_2D = np.nan*np.ones((maxGen-minGen+1,10*12)); HcondNE2_2D = np.nan*np.ones((maxGen-minGen+1,10*12)) 

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

    n1 = TE_EN_gen.shape[0]; print('For this gen, n = ',n1)

    # For slope
    HcondE = np.concatenate((HcondE,HcondE_gen),axis=0)
    HcondN = np.concatenate((HcondN,HcondN_gen),axis=0)
    TE_EN = np.concatenate((TE_EN,TE_EN_gen),axis=0)
    TE_NE = np.concatenate((TE_NE,TE_NE_gen),axis=0)
    TE_star_EN = np.concatenate((TE_star_EN,TE_star_EN_gen),axis=0)
    TE_star_NE = np.concatenate((TE_star_NE,TE_star_NE_gen),axis=0)
    TE_star2_EN = np.concatenate((TE_star2_EN,TE_star2_EN_gen),axis=0)
    TE_star2_NE = np.concatenate((TE_star2_NE,TE_star2_NE_gen),axis=0)
    HcondEN = np.concatenate((HcondEN,HcondEN_gen),axis=0)
    HcondNE = np.concatenate((HcondNE,HcondNE_gen),axis=0)
    HcondEN2 = np.concatenate((HcondEN2,HcondEN2_gen),axis=0)
    HcondNE2 = np.concatenate((HcondNE2,HcondNE2_gen),axis=0)

    # For box plots
    HcondE_2D[iGen-minGen,:len(HcondE_gen)] = HcondE_gen
    HcondN_2D[iGen-minGen,:len(HcondN_gen)] = HcondN_gen
    TE_EN_2D[iGen-minGen,:len(TE_EN_gen)] = TE_EN_gen
    TE_NE_2D[iGen-minGen,:len(TE_NE_gen)] = TE_NE_gen
    TE_star_EN_2D[iGen-minGen,:len(TE_star_EN_gen)] = TE_star_EN_gen
    TE_star_NE_2D[iGen-minGen,:len(TE_star_NE_gen)] = TE_star_NE_gen
    TE_star2_EN_2D[iGen-minGen,:len(TE_star2_EN_gen)] = TE_star2_EN_gen
    TE_star2_NE_2D[iGen-minGen,:len(TE_star2_NE_gen)] = TE_star2_NE_gen
    HcondEN_2D[iGen-minGen,:len(HcondEN_gen)] = HcondEN_gen
    HcondNE_2D[iGen-minGen,:len(HcondNE_gen)] = HcondNE_gen
    HcondEN2_2D[iGen-minGen,:len(HcondEN2_gen)] = HcondEN2_gen
    HcondNE2_2D[iGen-minGen,:len(HcondNE2_gen)] = HcondNE2_gen

    Genarr = np.concatenate((Genarr,iGen*np.ones(TE_EN_gen.shape)),axis=0)

print('Size of entire array = ',TE_EN.shape)

netTE_EN = 100*(TE_EN - TE_NE)
netnormTE_EN = 100*(TE_EN/HcondN - TE_NE/HcondE)
netTE_star_EN = 100*(TE_star_EN - TE_star_NE)
netnormTE_star_EN = 100*(TE_star_EN/HcondNE - TE_star_NE/HcondEN)
netTE_star2_EN = 100*(TE_star2_EN - TE_star2_NE)
netnormTE_star2_EN = 100*(TE_star_EN/HcondNE2 - TE_star_NE/HcondEN2)

x = Genarr.reshape(-1)
print('net TE ----> ')
netTE_TS = TS_slope(x,netTE_EN.reshape(-1))
print('net TE_C1 ----> ')
netTE_star_TS = TS_slope(x,netTE_star_EN.reshape(-1))
print('net TE_C2 ----> ')
netTE_star2_TS = TS_slope(x,netTE_star2_EN.reshape(-1))
print('net norm TE ----> ')
netnormTE_TS = TS_slope(x,netnormTE_EN.reshape(-1))
print('net norm TE_C1 ----> ')
netnormTE_star_TS = TS_slope(x,netnormTE_star_EN.reshape(-1))
print('net norm TE_C2 ----> ')
netnormTE_star2_TS = TS_slope(x,netnormTE_star2_EN.reshape(-1))

netTE_EN_2D = 100*(TE_EN_2D - TE_NE_2D)
netnormTE_EN_2D = 100*(TE_EN_2D/HcondN_2D - TE_NE_2D/HcondE_2D)
netTE_star_EN_2D = 100*(TE_star_EN_2D - TE_star_NE_2D)
netnormTE_star_EN_2D = 100*(TE_star_EN_2D/HcondNE_2D - TE_star_NE_2D/HcondEN_2D)
netTE_star2_EN_2D = 100*(TE_star2_EN_2D - TE_star2_NE_2D)
netnormTE_star2_EN_2D = 100*(TE_star2_EN_2D/HcondNE2_2D - TE_star2_NE_2D/HcondEN2_2D)

## -------

# Plot box-plots & mean for net TE, TE*, TE** of diff generations
data1 = netTE_EN_2D[0,:]; data1 = data1[~np.isnan(data1)]
data2 = netTE_EN_2D[1,:]; data2 = data2[~np.isnan(data2)]
data3 = netTE_EN_2D[2,:]; data3 = data3[~np.isnan(data3)]
data4 = netTE_EN_2D[3,:]; data4 = data4[~np.isnan(data4)]
data = [data1,data2,data3,data4] 

fig, (ax)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
c = "paleturquoise"; medc = "royalblue"
plt.boxplot(data, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c,linewidth=lw), #, color=c),
            capprops=dict(linewidth=lw),
            whiskerprops=dict(linewidth=lw),
            flierprops=dict(linewidth=lw),
            medianprops=dict(color=medc,linewidth=lw),
            )
plt.plot([0,5], [0,0], color='k', linestyle='dashed', linewidth = 1.0)

m = netTE_TS[0]; c = netTE_TS[1] # slope, intercept
x = np.array([2,3,4,5]); y = m*x+c
x2 = np.array([1.5,2,3,4,5,5.5]); y2 = m*x2+c
plt.plot(x2-1,y2,'-r')

fig.text(0.17, 0.83, r'TS $slope =$'+r"${0:.3g}$".format(m)+ \
            r', $p = $'+"{0:.2g}".format(netTE_TS[2]), color = "darkred")
ax.set_xticklabels(['2','3','4','5']); ax.set_xlabel('Generation')
ax.set_ylabel(r'$\mathrm{net \; TE}_{E → N} \; (\%)$')
ax.set_ylim([-4,5]); ax.set_yticks(np.arange(-4,6,1))
plt.savefig(fdir + 'boxplot_net_TE_EN_perc' + fstr, dpi=150, bbox_inches="tight")

## -------

data1B = netTE_star_EN_2D[0,:]; data1B = data1B[~np.isnan(data1B)]
data2B = netTE_star_EN_2D[1,:]; data2B = data2B[~np.isnan(data2B)]
data3B = netTE_star_EN_2D[2,:]; data3B = data3B[~np.isnan(data3B)]
data4B = netTE_star_EN_2D[3,:]; data4B = data4B[~np.isnan(data4B)]
dataB = [data1B,data2B,data3B,data4B] 

fig, (ax)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
c = "paleturquoise"; medc = "royalblue"
plt.boxplot(dataB, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c,linewidth=lw), #, color=c),
            capprops=dict(linewidth=lw),
            whiskerprops=dict(linewidth=lw),
            flierprops=dict(linewidth=lw),
            medianprops=dict(color=medc,linewidth=lw),
            )
plt.plot([0,5], [0,0], color='k', linestyle='dashed', linewidth = 1.0)

m = netTE_star_TS[0]; c = netTE_star_TS[1] # slope, intercept
x = np.array([2,3,4,5]); y = m*x+c
x2 = np.array([1.5,2,3,4,5,5.5]); y2 = m*x2+c
plt.plot(x2-1,y2,'-r')


fig.text(0.17, 0.83, r'TS $slope =$'+r"${0:.3g}$".format(m)+ \
            r', $p = $'+"{0:.2g}".format(netTE_star_TS[2]), color = "darkred")
ax.set_xticklabels(['2','3','4','5']); ax.set_xlabel('Generation')
ax.set_ylabel(r'$\mathrm{net \; TE}^{C}_{E → N} \; (\%)$')
ax.set_ylim([-4,5]); ax.set_yticks(np.arange(-4,6,1))
plt.savefig(fdir + 'boxplot_net_TE_C1_EN_perc' + fstr, dpi=150, bbox_inches="tight")

## -------

data1B = netTE_star2_EN_2D[0,:]; data1B = data1B[~np.isnan(data1B)]
data2B = netTE_star2_EN_2D[1,:]; data2B = data2B[~np.isnan(data2B)]
data3B = netTE_star2_EN_2D[2,:]; data3B = data3B[~np.isnan(data3B)]
data4B = netTE_star2_EN_2D[3,:]; data4B = data4B[~np.isnan(data4B)]
dataB = [data1B,data2B,data3B,data4B] 

fig, (ax)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
c = "paleturquoise"; medc = "royalblue"
plt.boxplot(dataB, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c,linewidth=lw), #, color=c),
            capprops=dict(linewidth=lw),
            whiskerprops=dict(linewidth=lw),
            flierprops=dict(linewidth=lw),
            medianprops=dict(color=medc,linewidth=lw),
            )
plt.plot([0,5], [0,0], color='k', linestyle='dashed', linewidth = 1.0)

m = netTE_star2_TS[0]; c = netTE_star2_TS[1] # slope, intercept
x = np.array([2,3,4,5]); y = m*x+c
x2 = np.array([1.5,2,3,4,5,5.5]); y2 = m*x2+c
plt.plot(x2-1,y2,'-r')

fig.text(0.17, 0.83, r'TS $slope =$'+r"${0:.3g}$".format(m)+ \
            r', $p = $'+"{0:.2g}".format(netTE_star2_TS[2]), color = "darkred")
ax.set_xticklabels(['2','3','4','5']); ax.set_xlabel('Generation')
ax.set_ylabel(r'$\mathrm{net \; TE}^{FC}_{E → N} \; (\%)$')
ax.set_ylim([-4,5]); ax.set_yticks(np.arange(-4,6,1))
plt.savefig(fdir + 'boxplot_net_TE_C2_EN_perc' + fstr, dpi=150, bbox_inches="tight")


## -------

# Plot box-plots for net norm TE of diff generations
data1 = netnormTE_EN_2D[0,:]; data1 = data1[~np.isnan(data1)]
data2 = netnormTE_EN_2D[1,:]; data2 = data2[~np.isnan(data2)]
data3 = netnormTE_EN_2D[2,:]; data3 = data3[~np.isnan(data3)]
data4 = netnormTE_EN_2D[3,:]; data4 = data4[~np.isnan(data4)]
data = [data1,data2,data3,data4] 

fig, (ax)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
c = "paleturquoise"; medc = "royalblue"
plt.boxplot(data, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c,linewidth=lw), #, color=c),
            capprops=dict(linewidth=lw),
            whiskerprops=dict(linewidth=lw),
            flierprops=dict(linewidth=lw),
            medianprops=dict(color=medc,linewidth=lw),
            )
plt.plot([0,5], [0,0], color='k', linestyle='dashed', linewidth = 1.0)

m = netnormTE_TS[0]; c = netnormTE_TS[1] # slope, intercept
x = np.array([2,3,4,5]); y = m*x+c
x2 = np.array([1.5,2,3,4,5,5.5]); y2 = m*x2+c
plt.plot(x2-1,y2,'-r')

fig.text(0.17, 0.83, r'TS $slope =$'+r"${0:.3g}$".format(m)+ \
            r', $p = $'+"{0:.2g}".format(netnormTE_TS[2]), color = "darkred")
ax.set_xticklabels(['2','3','4','5']); ax.set_xlabel('Generation')
ax.set_ylabel(r"$\mathrm{net \; TE^{'}}_{E → N} \; (\%)$")
ax.set_ylim([-4,5]); ax.set_yticks(np.arange(-4,6,1))
plt.savefig(fdir + 'boxplot_net_norm_TE_EN_perc' + fstr, dpi=150, bbox_inches="tight")


## -------

data1B = netnormTE_star_EN_2D[0,:]; data1B = data1B[~np.isnan(data1B)]
data2B = netnormTE_star_EN_2D[1,:]; data2B = data2B[~np.isnan(data2B)]
data3B = netnormTE_star_EN_2D[2,:]; data3B = data3B[~np.isnan(data3B)]
data4B = netnormTE_star_EN_2D[3,:]; data4B = data4B[~np.isnan(data4B)]
dataB = [data1B,data2B,data3B,data4B] 

fig, (ax)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
c = "paleturquoise"; medc = "royalblue"
plt.boxplot(dataB, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c,linewidth=lw), #, color=c),
            capprops=dict(linewidth=lw),
            whiskerprops=dict(linewidth=lw),
            flierprops=dict(linewidth=lw),
            medianprops=dict(color=medc,linewidth=lw),
            )
plt.plot([0,5], [0,0], color='k', linestyle='dashed', linewidth = 1.0)

m = netnormTE_star_TS[0]; c = netnormTE_star_TS[1] # slope, intercept
x = np.array([2,3,4,5]); y = m*x+c
x2 = np.array([1.5,2,3,4,5,5.5]); y2 = m*x2+c
plt.plot(x2-1,y2,'-r')

fig.text(0.17, 0.83, r'TS $slope =$'+r"${0:.3g}$".format(m)+ \
            r', $p = $'+"{0:.2g}".format(netnormTE_star_TS[2]), color = "darkred")
ax.set_xticklabels(['2','3','4','5']); ax.set_xlabel('Generation')
ax.set_ylabel(r"$\mathrm{net \; TE^{C1'}}_{E → N} \; (\%)$")
ax.set_ylim([-4,5]); ax.set_yticks(np.arange(-4,6,1))
plt.savefig(fdir + 'boxplot_net_norm_TE_C1_EN_perc' + fstr, dpi=150, bbox_inches="tight")


## -------

data1B = netnormTE_star2_EN_2D[0,:]; data1B = data1B[~np.isnan(data1B)]
data2B = netnormTE_star2_EN_2D[1,:]; data2B = data2B[~np.isnan(data2B)]
data3B = netnormTE_star2_EN_2D[2,:]; data3B = data3B[~np.isnan(data3B)]
data4B = netnormTE_star2_EN_2D[3,:]; data4B = data4B[~np.isnan(data4B)]
dataB = [data1B,data2B,data3B,data4B] 

fig, (ax)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
c = "paleturquoise"; medc = "royalblue"
plt.boxplot(dataB, notch=True, patch_artist=True,
            boxprops=dict(facecolor=c,linewidth=lw), #, color=c),
            capprops=dict(linewidth=lw),
            whiskerprops=dict(linewidth=lw),
            flierprops=dict(linewidth=lw),
            medianprops=dict(color=medc,linewidth=lw),
            )
plt.plot([0,5], [0,0], color='k', linestyle='dashed', linewidth = 1.0)

m = netnormTE_star2_TS[0]; c = netnormTE_star2_TS[1] # slope, intercept
x = np.array([2,3,4,5]); y = m*x+c
x2 = np.array([1.5,2,3,4,5,5.5]); y2 = m*x2+c
plt.plot(x2-1,y2,'-r')

fig.text(0.17, 0.83, r'TS $slope =$'+r"${0:.3g}$".format(m)+ \
            r', $p = $'+r"{0:.2g}".format(netnormTE_star2_TS[2]), color = "darkred")
ax.set_xticklabels(['2','3','4','5']); ax.set_xlabel('Generation')
ax.set_ylabel(r"$\mathrm{net \; TE^{C2'}}_{E → N} \; (\%)$")
ax.set_ylim([-4,5]); ax.set_yticks(np.arange(-4,6,1))
plt.savefig(fdir + 'boxplot_net_norm_TE_C2_EN_perc' + fstr, dpi=150, bbox_inches="tight")

# plt.show()    