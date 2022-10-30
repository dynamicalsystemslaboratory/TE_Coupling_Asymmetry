# This func calculates H(X_k|X_k-1) = I( X_k : Y_k-delta | X_k-1 )
# --> entropy of X given its past X
# m = embedding dimension (2,3,4)
import numpy as np
from compute_entropy import *
import time
    
def cond_entropy_sym_k(Xall,m,k): #,delta
    
    Xall = np.reshape(Xall,(len(Xall))) # make sure X, Y are (N,) arrays (not Nx1)
    
    ## m = 1 case : already symbolized / no need to symbolize
    if (m == 1):
        piXall = Xall
    ## m = 2 case
    elif (m == 2):
        piXall = np.diff(Xall) > 0
    ## m = 3,4 cases
    else:
        piXall = symbolize_data(Xall,m)
    
    # Remove Nans from both time series
    piXall = piXall[np.isnan(piXall)==False]
    
    # X_n-1 = X, X_n(k) = X1, X_n+1 = X2
    N = len(piXall)     

    piXall = piXall.reshape((N,1)) # make sure X, Y are (N,1) arrays (not (N,))
                                   # required for np.append ...
    nst = k

    piX1 = np.array([]); piX2 = np.array([])
    for i in range(1,k+1):
        if (~np.any(piX1)): 
            piX1 = piXall[nst-i+1:N-i]
        else:
            piX1 = np.append(piX1,piXall[nst-i+1:N-i],axis=1)
            
    piX2 = piXall[nst+1:N]

    HX1X2 = compute_entropy(np.append(piX1,piX2,axis=1))
    HX1 = compute_entropy(piX1) # make sure each column is a timeseries

    # Compute cond entropy
    HcondX = HX1X2 - HX1

    return HcondX

def symbolize_data(X,m): 
    if (m == 3):
        # 3-histories
        X0 = X[0:X.shape[0]-2] # X_t
        X1 = X[1:X.shape[0]-1] # X_t+1
        X2 = X[2:X.shape[0]] # X_t+2
        symX = np.argsort( np.transpose(np.array([X0,X1,X2])), axis=1 ) # make sure each time series is a column, then sort along each row
        l,piX = np.unique(symX,axis=0,return_inverse=True)

    elif (m == 4):
        # 4-histories
        X0 = X[0:X.shape[0]-3] # X_t
        X1 = X[1:X.shape[0]-2] # X_t+1
        X2 = X[2:X.shape[0]-1] # X_t+2
        X3 = X[3:X.shape[0]] # X_t+3
        symX = np.argsort( np.transpose(np.array([X0,X1,X2,X3])), axis=1 ) # make sure each time series is a column, then sort along each row
        l,piX = np.unique(symX,axis=0,return_inverse=True)
    else:
        print('Sorry, we have not coded for m>4 yet! ')
    
    return piX