import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import basinhopping
from scipy.stats import norm
import pandas as pd

from transfer_entropy2 import *

# parameters
N = 1 # no. of simulations 
Nt = 10000000 # no. of time steps in each simulation
sampN = 1 # skip these many, i.e. downsampling fac 
sigma = 0.1 #0.1 
sigw = sigma; sigv = sigma

a = 0.3    # 0.3
b = 0.4    # 0.2
eps = 0.1  # 0.1
# alph = 2.0 # 2.0

fdir = './figures_bin/'; 
fdir2 = './files_bin2/'; 
fstr = 'a0_3_b0_4_eps0_1_sig0_1_r0_5_k1_l1_Nt10mil'

k = 1; l = 1
r = 0.5

alpharr = np.arange(0,10.5,0.5)
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

with open(fdir + 'readme_' + fstr + '.txt', 'w') as f:
    f.write('alph     # x or y out of bounds             x-range             y-range  \n')


Tx_y = np.zeros((Nalph,N)); Ty_x = np.zeros((Nalph,N))
Tx_y_star = np.zeros((Nalph,N)); Ty_x_star = np.zeros((Nalph,N))

print(" k = ", k, ".........")

for ialph in range(0,Nalph):
    alph = alpharr[ialph]
    print(" alpha = ", alph, "---------------------")

    for ksim in range(0,N):
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
            
        setzero = np.sum(x>1) + np.sum(x<0) + np.sum(y<0) + np.sum(y>1)
        print('# times x or y was set to 0 = '+ str(setzero) )

        with open(fdir + 'readme_' + fstr + '.txt', 'a') as f:
            f.write( str(alph) + '           ' + str(setzero) + '                     ' + str(np.float32(np.min(x))) + ': ' + str(np.float32(np.max(x))) + '  ' + str(np.float32(np.min(y))) + ': ' + str(np.float32(np.max(y))) + '\n')

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

        # Calc all
        Ty_x_all = transfer_entropy2(x,y,xmax,xmin,r,ymax,ymin,r,k,l)
        Tx_y_all = transfer_entropy2(y,x,xmax,xmin,r,ymax,ymin,r,k,l)

        Ty_x[ialph,ksim] = Ty_x_all[0]
        Tx_y[ialph,ksim] = Tx_y_all[0]

        Ty_x_star[ialph,ksim] = Ty_x_all[1]
        Tx_y_star[ialph,ksim] = Tx_y_all[1]

        print(Tx_y_all-Ty_x_all)

print('Loops done! Avg over sims...')

# Mean
meanTE_xy = np.mean(Tx_y,axis=1)
meanTE_star_xy = np.mean(Tx_y_star,axis=1) 
meanTE_yx = np.mean(Ty_x,axis=1)
meanTE_star_yx = np.mean(Ty_x_star,axis=1)

# Datasave
print('Datasave ....')
df = pd.DataFrame()
df['alpha']=alpharr
df['TExy']=meanTE_xy
df['TEyx']=meanTE_yx
df['TECxy']=meanTE_star_xy
df['TECyx']=meanTE_star_yx
df.to_csv(fdir2 + 'TE_vs_alph_' + fstr + '.csv')
print('Plot ....')

## Plot cond entropy for using X as target
plt.figure()
plt.plot(alpharr,meanTE_xy-meanTE_xy, linewidth=2, marker='*', markersize=10, label='netTE_{x → y}')
plt.plot(alpharr,meanTE_star_xy-meanTE_star_xy, linewidth=2, marker='*', markersize=10, label='netTE^{C}_{x → y}')
plt.title('a = '+str(a)+', b = '+str(b)+', eps = '+str(eps)+', r = '+str(r)+', k = '+str(k)+', l = '+str(l))
plt.xlabel(r'$\alpha$',fontsize=15)
plt.ylabel(r'$\mathrm{net \; TE}_{x → y}$',fontsize=15)
plt.legend(fontsize=12)
plt.savefig(fdir + 'netTE_all_vs_alph_' + fstr + '.png')

print('All done!')
