import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import basinhopping
from scipy.stats import norm
import pandas as pd

# from transfer_entropy_bin_all import *
from transfer_entropy3 import *

# parameters
simN = 3
N = 30 # no. of simulations 
Nt = 10000000 # no. of time steps in each simulation
sampN = 1 # skip these many, i.e. downsampling fac 
sigma = 0.1 #0.1 
sigw = sigma; sigv = sigma

# a = 0.4    # 0.3
# b = 0.5    # 0.2
eps = 0.1  # 0.1
# alph = 2.0 # 2.0

k = 1; l = 1
r = 0.5; rstr = '0_5'

# fdir = './figures_bin/'; 
fdir2 = './files_bin_sims2/'; 
fstr = 'eps0_1_sig0_1_r' + rstr + '_Nt10mil_'

alpharr = np.array([0,0.5,1,2,4,6,8])
Nalph = np.size(alpharr)

# Tent map 
mu=2.0
def tent(xx,mu):   
    F=0
    if ((xx > 0 ) & (xx <= 0.5)):
        F = mu*xx
    elif ((xx > 0.5 ) & (xx <= 1.0)):
        F = mu*(1-xx)
    return F

## Main loop ======================================================================
Tx_y = np.zeros((Nalph,N)); Ty_x = np.zeros((Nalph,N))
Tx_y_star = np.zeros((Nalph,N)); Ty_x_star = np.zeros((Nalph,N))
Tx_y_star2 = np.zeros((Nalph,N)); Ty_x_star2 = np.zeros((Nalph,N))

# Randomly pick a,b
aarr = np.random.uniform(0,0.45,N) # Later maybe try (-0.45, 0.45)
barr = np.random.uniform(0,0.45,N)

for ksim in range(0,N):
    a = aarr[ksim]; b = barr[ksim]
    print(" a = ", a, ", b = " , b , "---------------------")
    print(ksim)
        
    for ialph in range(0,Nalph):
        alph = alpharr[ialph]
        print(" alpha = ", alph, "---------------------")
        # Noise
        w = np.random.randn(Nt)
        v = np.random.randn(Nt)
        # initial condition
        x = np.zeros((Nt))
        y = np.zeros((Nt))
        x[0] = np.random.uniform(0,1)
        y[0] = np.random.uniform(0,1)

        # Nonlinear - Chaotic tent map
        for n in range(0,Nt-1):
            # print(x[n]); time.sleep(0.1)
            x[n+1] = a*tent(x[n],mu) + b*tent(y[n],mu) + sigw*w[n]
            y[n+1] = (b+eps)*tent(x[n],mu) + (a+(alph*eps))*tent(y[n],mu) + sigv*v[n]
            
        # In case solution blows up
        if (  ( (np.sum(np.isnan(x)*1==1)>0) or (np.sum(np.isnan(y)*1==1)>0) ) or (np.abs(np.max([x,y]))>1000) ):
            Ty_x[ialph,ksim] = -1000.0; Tx_y[ialph,ksim] = -1000.0
            Ty_x_star[ialph,ksim] = -1000.0; Tx_y_star[ialph,ksim] = -1000.0
            Ty_x_star2[ialph,ksim] = -1000.0; Tx_y_star2[ialph,ksim] = -1000.0
            continue

        setzero = np.sum(x>1) + np.sum(x<0) + np.sum(y<0) + np.sum(y>1)
        print('# times x or y was set to 0 = '+ str(setzero) )

        # Donwsample
        y = y[0:Nt:sampN]
        x = x[0:Nt:sampN]
        print('x : ' + str(np.min(x)) + ': ' + str(np.max(x)) + ', y :' + str(np.float32(np.min(y))) + ': ' + str(np.float32(np.max(y))) )

        #Normalize data
        x = (x-np.mean(x))/np.std(x)
        y = (y-np.mean(y))/np.std(y)
        print('x_norm : ' + str(np.min(x)) + ': ' + str(np.max(x)) + ', y_norm :' + str(np.float32(np.min(y))) + ': ' + str(np.float32(np.max(y))) )

        xmax = (np.max(x)); xmin = (np.min(x))
        ymax = (np.max(y)); ymin = (np.min(y))

        # TE - all
        Ty_x_all = transfer_entropy3(x,y,xmax,xmin,r,ymax,ymin,r,k,l)
        Tx_y_all = transfer_entropy3(y,x,xmax,xmin,r,ymax,ymin,r,k,l)

        # TE 
        Ty_x[ialph,ksim] = Ty_x_all[0]
        Tx_y[ialph,ksim] = Tx_y_all[0]

        # TE_C1
        Ty_x_star[ialph,ksim] = Ty_x_all[1]
        Tx_y_star[ialph,ksim] = Tx_y_all[1]

        # TE_C2
        Ty_x_star2[ialph,ksim] = Ty_x_all[2]
        Tx_y_star2[ialph,ksim] = Tx_y_all[2]

        print(Tx_y_all-Ty_x_all)

    # Datasave
    df = pd.DataFrame()
    df['alpha']=alpharr
    df['TExy']=Tx_y[:,ksim]
    df['TEyx']=Ty_x[:,ksim]
    df['TECxy']=Tx_y_star[:,ksim]
    df['TECyx']=Ty_x_star[:,ksim]
    df['TEC2xy']=Tx_y_star2[:,ksim]
    df['TEC2yx']=Ty_x_star2[:,ksim]
    df.to_csv(fdir2 + 'TE_vs_alph_' + fstr + str((simN-1)*N + ksim+1) + '.csv')

print('All done!')
