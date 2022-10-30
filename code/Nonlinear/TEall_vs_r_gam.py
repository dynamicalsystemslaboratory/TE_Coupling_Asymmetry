import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import basinhopping
from scipy.stats import norm
import pandas as pd

from transfer_entropy3 import *

# parameters
N = 1 # no. of simulations -- code is written only for N=1 !!!!!!!
Nt = 10000000 # no. of time steps in each simulation
sampN = 1 # skip these many, i.e. downsampling fac 
sigma = 0.1 #0.1 
sigw = sigma; sigv = sigma

a = 0.15    # 0.3 #### (0.1,0.2), (0.3,0.4), (0.4,0.3),  (0.5,0.1), (0.1,0.5)
b = 0.2    # 0.4
eps = 0.1  # 0.1
alph = 2.0 # 2.0

fdir = './figures_bin/'; 
fdir2 = './files_bin/'; 
fstr = 'a0_15_b0_2_eps0_1_alph2_sig0_1_Nt10mil'

rarr = np.array([0.05, 0.1, 0.2, 0.5]) # Norm
Nr = np.size(rarr)

k = 1; l = 1
gamarr = np.array([0.0,0.01,0.05,0.1,0.2,0.5]) 
Ngam = np.size(gamarr)
gamstr = np.array(['0','0_01','0_05','0_1','0_2','0_5'])

## Main loop ======================================================================
X = np.zeros((Nt,N)); Y = np.zeros((Nt,N))

for igam in range(0,Ngam):
    gam = gamarr[igam]
    print(" gamma = ", gam, "---------------------")

    Tx_y = np.zeros((Nr,N)); Ty_x = np.zeros((Nr,N))
    Tx_y_star = np.zeros((Nr,N)); Ty_x_star = np.zeros((Nr,N))
    Tx_y_star2 = np.zeros((Nr,N)); Ty_x_star2 = np.zeros((Nr,N))

    for ksim in range(0,N):
        # print(ksim)
        # Noise
        w = np.random.randn(Nt)
        v = np.random.randn(Nt)
        # initial condition
        x = np.zeros((Nt))
        y = np.zeros((Nt))
        x[0] = np.random.randn(1)
        y[0] = np.random.randn(1)

        # Non-linear
        for n in range(0,Nt-1):
            x[n+1] = a*x[n] + b*y[n] + gam*(x[n]**3) + sigw*w[n]
            y[n+1] = (b+eps)*x[n] + (a+(alph*eps))*y[n] + gam*(y[n]**3) + sigv*v[n]
            
        # Donwsample
        y = y[0:Nt:sampN]
        x = x[0:Nt:sampN]
        print('x : ' + str(np.min(x)) + ': ' + str(np.max(x)) + ', y :' + str(np.float32(np.min(y))) + ': ' + str(np.float32(np.max(y))) )

        #Normalize data
        x = (x-np.mean(x))/np.std(x)
        y = (y-np.mean(y))/np.std(y)
        print('x_norm : ' + str(np.min(x)) + ': ' + str(np.max(x)) + ', y_norm :' + str(np.float32(np.min(y))) + ': ' + str(np.float32(np.max(y))) )

        X[:,ksim] = x; Y[:,ksim] = y; 

    for ir in range(0,Nr):
        binwidth = rarr[ir]
        print(" r = ", binwidth, "---------------------")

        for ksim in range(0,N):
            x = X[:,ksim]
            y = Y[:,ksim]

            xmax = (np.max(x)); xmin = (np.min(x))
            ymax = (np.max(y)); ymin = (np.min(y))

            # Calc all
            Ty_x_all = transfer_entropy3(x,y,xmax,xmin,binwidth,ymax,ymin,binwidth,k,l)
            Tx_y_all = transfer_entropy3(y,x,xmax,xmin,binwidth,ymax,ymin,binwidth,k,l)

            Ty_x[ir,ksim] = Ty_x_all[0]
            Tx_y[ir,ksim] = Tx_y_all[0]

            Ty_x_star[ir,ksim] = Ty_x_all[1]
            Tx_y_star[ir,ksim] = Tx_y_all[1]
            
            Ty_x_star2[ir,ksim] = Ty_x_all[2]
            Tx_y_star2[ir,ksim] = Tx_y_all[2]

            print(Tx_y_all-Ty_x_all)

    # Mean
    meanTE_xy = np.mean(Tx_y,axis=1)
    meanTE_star_xy = np.mean(Tx_y_star,axis=1)
    meanTE_star2_xy = np.mean(Tx_y_star2,axis=1) 
    meanTE_yx = np.mean(Ty_x,axis=1)
    meanTE_star_yx = np.mean(Ty_x_star,axis=1)
    meanTE_star2_yx = np.mean(Ty_x_star2,axis=1)

    # Datasave
    print('Datasave ....')
    fstr2= '_gam' + gamstr[igam]
    df = pd.DataFrame()
    df['r']=rarr
    df['TExy']=meanTE_xy
    df['TEyx']=meanTE_yx
    df['TECxy']=meanTE_star_xy
    df['TECyx']=meanTE_star_yx
    df['TEC2xy']=meanTE_star2_xy
    df['TEC2yx']=meanTE_star2_yx
    df.to_csv(fdir2 + 'TE_bin_all_vs_r_' + fstr + fstr2 + '.csv')

print('All done!')
