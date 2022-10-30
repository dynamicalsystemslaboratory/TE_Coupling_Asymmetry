import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd

from transfer_entropy_kde_all import *
from numpy import loadtxt

plt.rcParams.update({
    "font.size" : 17,
    "lines.linewidth" : 2.0,
    "font.family": "arial"
    })

data = loadtxt("./SantaFe_data/b1.txt",usecols = (0,1,2))

# Data is at a sampling frequency of 0.5 s
hr1 = data[:,0] # Heart-rate
resp1 = data[:,1] # Respiration force (same as breath rate)
dt = 0.5

data = loadtxt("./SantaFe_data/b2.txt",usecols = (0,1,2))

# Data is at a sampling frequency of 0.5 s
hr2 = data[:,0] # Heart-rate
resp2 = data[:,1] # Respiration force (same as breath rate)

hr = np.concatenate((hr1,hr2),axis=0)
resp = np.concatenate((resp1,resp2),axis=0)

# Schreiber used data samples 2350–3550 only
Nstart = 0 #2350
Nend = 34000+1 #3550+1
x = hr[Nstart:Nend]
y = resp[Nstart:Nend]

Nt= (Nend-1)-Nstart;  t0=0.0
tarr = t0 + dt*np.arange(Nt)

#Normalize data
mux = np.mean(x)
muy = np.mean(y)
sigmax = np.std(x)
sigmay = np.std(y)
x = (x-mux)/sigmax
y = (y-muy)/sigmay
print('After normalization, x = ' +str(np.mean(x))+' +/- '+str(np.std(x))+' and y = '+str(np.mean(y))+' +/- '+str(np.std(y)))

# Plot data
fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=(7, 2.5))
ax1.plot(tarr[400:1800],y[400:1800], color = "b", label = "breath rate")
ax1.set_xlabel(r"$t$ (seconds)"); ax1.set_ylabel(r"breath rate ($BR$)")
ax1.set_ylim([-4,4])
ax1.set_xlim([200,900])
fig.savefig("./figures2/normBR_v_t.eps", dpi=150, bbox_inches="tight")

fig, (ax1)  = plt.subplots(1, 1,sharey='row',figsize=(7, 2.5))
ax1.plot(tarr[400:1800],x[400:1800], color = "k", label = "heart rate")
ax1.set_xlabel(r"$t$ (seconds)"); ax1.set_ylabel(r"heart rate ($HR$)")
ax1.set_ylim([-2,2])
ax1.set_xlim([200,900])
fig.savefig("./figures2/normHR_v_t.eps", dpi=150, bbox_inches="tight")

print("Paused!")
time.sleep(50000)

# # parameters
sampN = 1 # skip these many, i.e. downsampling fac 

rarr= np.array([0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,1.0])#,0.5,1.0])
Nr = np.size(rarr)

MIx_y = np.zeros(Nr); MIy_x = np.zeros(Nr)
Tx_y = np.zeros(Nr); Ty_x = np.zeros(Nr)
Tx_y_star = np.zeros(Nr); Ty_x_star = np.zeros(Nr)
Tx_y_star2 = np.zeros(Nr); Ty_x_star2 = np.zeros(Nr)

## Main loop ======================================================================

for ir in range(0,Nr):
    bandwidth = rarr[ir]
    print(" r = ", bandwidth, "---------------------")

    # Donwsample
    y = y[0:Nt:sampN]
    x = x[0:Nt:sampN]

    # Calc all
    Ty_x_all = transfer_entropy_kde_all(x,y,bandwidth,'tophat') 
    Tx_y_all = transfer_entropy_kde_all(y,x,bandwidth,'tophat') 

    MIy_x[ir] = Ty_x_all[0] # MI
    MIx_y[ir] = Tx_y_all[0]

    Ty_x[ir] = Ty_x_all[1] # TE
    Tx_y[ir] = Tx_y_all[1]

    Ty_x_star[ir] = Ty_x_all[2] # TE_C1
    Tx_y_star[ir] = Tx_y_all[2]
    
    Ty_x_star2[ir] = Ty_x_all[3] # TE_C2
    Tx_y_star2[ir] = Tx_y_all[3]

    print('net MI, TE, TE*, TE** = ' + str(MIx_y[ir]-MIy_x[ir]) + ', ' + str(Tx_y[ir]-Ty_x[ir]) + ', ' + str(Tx_y_star[ir]-Ty_x_star[ir]) + ', ' + str(Tx_y_star2[ir]-Ty_x_star2[ir]))


df = pd.DataFrame()
df['rarr']=rarr
df['MIx_y']=MIx_y; df['MIy_x']=MIy_x
df['Tx_y']=Tx_y; df['Ty_x']=Ty_x
df['Tx_y_star']=Tx_y_star; df['Ty_x_star']=Ty_x_star
df['Tx_y_star2']=Tx_y_star2; df['Ty_x_star2']=Ty_x_star2

# Modify string names for writing
fdir= './files/'; 
df.to_csv(fdir + 'TE_all_1.csv')