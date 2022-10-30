import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import matplotlib as mpl

plt.rcParams.update({
    "font.size" : 17,
    "lines.linewidth" : 2.0,
    "font.family": "arial"
    })

fsz = (7,3.5)

## Load data
fdir= './files/'; 
str = fdir + 'TE_all_1.csv'
df = pd.read_csv(str, index_col=0)
print(df.head())

rarr = np.array(df['rarr'])
MIx_y = np.array(df['MIx_y']); MIy_x = np.array(df['MIy_x'])
Tx_y = np.array(df['Tx_y']); Ty_x = np.array(df['Ty_x'])
Tx_y_star = np.array(df['Tx_y_star']); Ty_x_star = np.array(df['Ty_x_star'])
Tx_y_star2 = np.array(df['Tx_y_star2']); Ty_x_star2 = np.array(df['Ty_x_star2'])

## Plot TE, MI, TE*, TE** vs r  to compare X->Y and Y->X 

fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
ax1.plot(rarr,MIx_y,color='g', marker='*', markersize=10, label=r'$\mathrm{MI}_{HR → BR}$')
ax1.plot(rarr,MIy_x,color='g', linestyle='dashed',label=r'$\mathrm{MI}_{BR → HR}$')
ax1.set_xscale('log')
ax1.set_xlabel(r'$r$'); ax1.set_ylabel(r'$\mathrm{MI}$')
ax1.legend()
plt.savefig('./figures2/MI_vs_r_KDE.eps', dpi=150, bbox_inches="tight")

fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
ax1.plot(rarr,Tx_y,color='k', marker='*', markersize=10, label=r'$\mathrm{TE}_{HR → BR}$')
ax1.plot(rarr,Ty_x,color='k', linestyle='dashed',label=r'$\mathrm{TE}_{BR → HR}$')
ax1.set_xscale('log')
ax1.set_ylim([0,5])
ax1.set_xlabel(r'$r$'); ax1.set_ylabel(r'$\mathrm{TE}$')
ax1.legend()
plt.savefig('./figures2/TE_vs_r_KDE.eps', dpi=150, bbox_inches="tight")

fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
ax1.plot(rarr,Tx_y_star,color='r', marker='*', markersize=10, label=r'$\mathrm{TE}^{C1}_{HR → BR}$')
ax1.plot(rarr,Ty_x_star,color='r', linestyle='dashed',label=r'$\mathrm{TE}^{C1}_{BR → HR}$')
ax1.set_xscale('log')
ax1.set_ylim([0.0,2.0])
ax1.set_xlabel(r'$r$'); ax1.set_ylabel(r'$\mathrm{TE}^{C1}$')
ax1.legend()
plt.savefig('./figures2/TEstar_vs_r_KDE.eps', dpi=150, bbox_inches="tight")


fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=fsz)
ax1.plot(rarr,Tx_y_star2,color='b', marker='*', markersize=10, label=r'$\mathrm{TE}^{C2}_{HR → BR}$')
ax1.plot(rarr,Ty_x_star2,color='b', linestyle='dashed',label=r'$\mathrm{TE}^{C2}_{BR → HR}$')
ax1.set_xscale('log')
ax1.set_ylim([0.0,1.25])
ax1.set_xlabel(r'$r$'); ax1.set_ylabel(r'$\mathrm{TE}^{C2}$')
ax1.legend()
plt.savefig('./figures2/TEstar2_vs_r_KDE.eps', dpi=150, bbox_inches="tight")

fig=plt.figure()
fig.set_size_inches(7, 3.5)
plt.plot(rarr,MIx_y,color='g', linewidth=2, marker='*', label='MIx_y')
plt.plot(rarr,MIy_x,color='g', linewidth=2,linestyle='dashed',label='MIy_x')
plt.plot(rarr,Tx_y,color='k', linewidth=2, marker='*', label='TEx_y')
plt.plot(rarr,Ty_x,color='k', linewidth=2,linestyle='dashed',label='TEy_x')
plt.plot(rarr,Tx_y_star,color='r', linewidth=2, marker='*', label='TE*x_y')
plt.plot(rarr,Ty_x_star,color='r', linewidth=2,linestyle='dashed',label='TE*y_x')
plt.plot(rarr,Tx_y_star2,color='b', linewidth=2, marker='*', label='TE**x_y')
plt.plot(rarr,Ty_x_star2,color='b', linewidth=2,linestyle='dashed',label='TE**y_x')
plt.xscale('log')
plt.xlabel('r',fontsize=15)
plt.legend(fontsize=15)
plt.savefig('./figures/AllTE_vs_r_KDE.eps')

fig =plt.figure()
fig.set_size_inches(7, 3.5)
plt.plot(rarr,MIx_y,color='g', linewidth=2, marker='*', label='MIx_y')
plt.plot(rarr,MIy_x,color='g', linewidth=2,linestyle='dashed',label='MIy_x')
plt.plot(rarr,Tx_y,color='k', linewidth=2, marker='*', label='TEx_y')
plt.plot(rarr,Ty_x,color='k', linewidth=2,linestyle='dotted',label='TEy_x')
plt.xscale('log')
plt.xlim([0.01,1.0]); plt.ylim([0,5])
plt.legend(fontsize=15)
plt.savefig('./figures/MI_TE_vs_r_KDE.eps')
# plt.show()