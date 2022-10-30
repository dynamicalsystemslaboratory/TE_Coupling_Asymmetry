## Computes TE,TE_C using binning
import numpy as np
import pandas as pd
from compute_entropy import *
import math

def symbolize(X,min,max,nbins):
    binedges = np.linspace(min, max, num=nbins+1)
    symbols = np.linspace(0, nbins, num=nbins, endpoint=False,dtype=int)

    df= pd.DataFrame({'ts': X})
    df['syms'] = pd.cut(x=df['ts'], bins=binedges,labels=symbols)
    symarr = df['syms'].to_numpy()
    return symarr

def transfer_entropy2(Xall,Yall,xmax,xmin,xbinwidth,ymax,ymin,ybinwidth,k,l):  
    # This func calculates 
    # TE_Y_to_X = I( X_n+1 : Y_n(l) | X_n(k) )
    # TE_C_Y_to_X = I( X_n+1 : Y_n(l) | X_n(k) , Y_n-l(l) ) ...

    # xnbins, ynbins: # of equally spaced points where we want to determine PDF to compute TE
    # k = # of time histories of target X to condition on
    # l = # of time histories of source Y to condition on

    TY_X = np.zeros(2)  # TE, TE_C

    Xall = np.reshape(Xall,(len(Xall))) # make sure X, Y are (N,) arrays (not Nx1) 
    Yall = np.reshape(Yall,(len(Yall))) # make sure X, Y are (N,) arrays (not Nx1) 
    N = len(Xall)

    # Symbolize time series by binning
    xnbins = int(math.ceil((xmax-xmin)/xbinwidth))
    ynbins = int(math.ceil((ymax-ymin)/ybinwidth))
    x_adjust=(xnbins*xbinwidth)-(xmax-xmin) # adjust xmin and max to keep binsize fixed
    xmax = xmax + (x_adjust/2.0); xmin = xmin - (x_adjust/2.0)
    y_adjust=(ynbins*ybinwidth)-(ymax-ymin)
    ymax = ymax + (y_adjust/2.0); ymin = ymin - (y_adjust/2.0)

    piXall = symbolize(Xall,xmin,xmax,xnbins)
    piYall = symbolize(Yall,ymin,ymax,ynbins)

    piXall = piXall.reshape((N,1)) # make sure X, Y are (N,1) arrays (not (N,))
    piYall = piYall.reshape((N,1)) 

    # X_n(2k) = XX1, X_n(k) = X1, X_n+1 = X2
    # Y_n(2l) = YY1, Y_n-l(l) = Y, Y_n(l) = Y1
    nst = np.max([2*l-1,2*k-1])

    piXX1 = np.array([]); piX1 = np.array([]); piX2 = np.array([])
    piYY1 = np.array([]); piY = np.array([]); piY1 = np.array([])

    for i in range(1,2*k+1):
        if (~np.any(piXX1)): # If empty create one
            piXX1 = piXall[nst-i+1:N-i]
        else:
            piXX1 = np.append(piXX1,piXall[nst-i+1:N-i],axis=1)

    for i in range(1,k+1):
        if (~np.any(piX1)): 
            piX1 = piXall[nst-i+1:N-i]
        else:
            piX1 = np.append(piX1,piXall[nst-i+1:N-i],axis=1)

    piX2 = piXall[nst+1:N]

    for i in range(1,2*l+1):
        if (~np.any(piYY1)): 
            piYY1 = piYall[nst-i+1:N-i]
        else:
            piYY1 = np.append(piYY1,piYall[nst-i+1:N-i],axis=1)

    for i in range(1,l+1):
        if (~np.any(piY)): 
            piY = piYall[nst-i+1-l:N-i-l]
        else:
            piY = np.append(piY,piYall[nst-i+1-l:N-i-l],axis=1)
            
        if (~np.any(piY1)): 
            piY1 = piYall[nst-i+1:N-i]
        else:
            piY1 = np.append(piY1,piYall[nst-i+1:N-i],axis=1)

    # For TE_C
    HYX1Y1X2 = compute_entropy(np.concatenate((piYY1,piX1,piX2),axis=1))
    HYX1Y1 = compute_entropy(np.concatenate((piYY1,piX1),axis=1))
    HYX1X2 = compute_entropy(np.concatenate((piY,piX1,piX2),axis=1))
    HYX1 = compute_entropy(np.concatenate((piY,piX1),axis=1))

    # For TE
    HX1Y1X2 = compute_entropy(np.concatenate((piX1,piY1,piX2),axis=1))
    HX1Y1 = compute_entropy(np.concatenate((piX1,piY1),axis=1))
    HX1X2 = compute_entropy(np.concatenate((piX1,piX2),axis=1))
    HX1 = compute_entropy(piX1)

    # Compute TE_C
    TY_X[1] = HYX1X2 - HYX1 - HYX1Y1X2 + HYX1Y1
    # Compute TE
    TY_X[0] = HX1X2 - HX1 - HX1Y1X2 + HX1Y1

    return TY_X