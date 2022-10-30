## TE based on symbolic approach
import numpy as np
from compute_entropy import *
# import time
    
def transfer_entropy_all_Hcond_sym_k(Xall,Yall,m,k,l):
    # This func calculates 
    # TE_Y_to_X = I( X_n+1 : Y_n | X_n(k) )
    # TE*_Y_to_X = I( X_n+1 : Y_n | X_n(k) , Y_n-1 )

    # m = embedding dimension (2,3,4)
    # k = # of time histories of X to condition on
    # l = # of time histories of Y to condition on

    TY_X = np.zeros(6)  # MI, TE, TE_C1, TE_C2, HcondXY, HcondXY2

    Xall = np.reshape(Xall,(len(Xall))) # make sure X, Y are (N,) arrays (not Nx1)  
    Yall = np.reshape(Yall,(len(Yall))) 
    
    ## m = 1 case : already symbolized / no need to symbolize
    if (m == 1):
        piXall = Xall
        piYall = Yall
    ## m = 2 case
    elif (m == 2):
        piXall = np.diff(Xall) > 0
        piYall = np.diff(Yall) > 0
    ## m = 3,4 cases
    else:
        piXall = symbolize_data(Xall,m)
        piYall = symbolize_data(Yall,m)
    
    # Remove entries from both t.s. if 
    # the value of any of the two t.s. is nan
    piXall = piXall[np.isnan(piXall)==False]
    piYall = piYall[np.isnan(piXall)==False]

    piYall = piYall[np.isnan(piYall)==False]
    piXall = piXall[np.isnan(piYall)==False]

    N = piXall.shape[0]   
    piXall = piXall.reshape((N,1)) # make sure X, Y are (N,1) arrays (not (N,))
    piYall = piYall.reshape((N,1)) # required for np.append ...
    
    # X_n(k+1) = XX1, X_n(k) = X1, X_n+1 = X2
    # Y_n(l+1) = YY1, Y_n-1(l) = Y, Y_n(l) = Y1
    nst = np.max([l,k])

    piXX1 = np.array([]);piX1 = np.array([]);piX2 = np.array([])
    piYY1 = np.array([]);piY = np.array([]);piY1 = np.array([])
    for i in range(1,k+1+1):
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

    for i in range(1,l+1+1):
        if (~np.any(piYY1)): 
            piYY1 = piYall[nst-i+1:N-i]
        else:
            piYY1 = np.append(piYY1,piYall[nst-i+1:N-i],axis=1)

    for i in range(1,l+1):
        if (~np.any(piY)): 
            piY = piYall[nst-i:N-1-i]
        else:
            piY = np.append(piY,piYall[nst-i:N-1-i],axis=1)
        if (~np.any(piY1)): 
            piY1 = piYall[nst-i+1:N-i]
        else:
            piY1 = np.append(piY1,piYall[nst-i+1:N-i],axis=1)

    # p(Y_(k-1), X_(k), Y_(k), X_(k+1))
    HYX1X2 = compute_entropy(np.append(np.append(piY,piX1,axis=1),piX2,axis=1))
    HYX1 = compute_entropy(np.append(piY,piX1,axis=1))
    HYY1X1X2 = compute_entropy(np.append(piYY1,np.append(piX1,piX2,axis=1),axis=1))
    HYY1X1 = compute_entropy(np.append(piYY1,piX1,axis=1)) 
    
    HXX1YX2 = compute_entropy(np.append(piXX1,np.append(piY,piX2,axis=1),axis=1))
    HXX1Y = compute_entropy(np.append(piXX1,piY,axis=1)) 
    HXX1YY1X2 = compute_entropy(np.append(piXX1,np.append(piYY1,piX2,axis=1),axis=1))
    HXX1YY1 = compute_entropy(np.append(piXX1,piYY1,axis=1)) 
    
    HX1X2 = compute_entropy(np.append(piX1,piX2,axis=1))
    HX1 = compute_entropy(piX1) # make sure each column is a timeseries
    HX1Y1X2 = compute_entropy(np.append(piX1,np.append(piY1,piX2,axis=1),axis=1))
    HX1Y1 = compute_entropy(np.append(piX1,piY1,axis=1)) 
    
    HX2 = compute_entropy(piX2)
    HY1X2 = compute_entropy(np.append(piY1,piX2,axis=1))
    HY1 = compute_entropy(piY1)

    # Compute TE_C2
    TY_X[3] = HXX1YX2 - HXX1Y - HXX1YY1X2 + HXX1YY1
    # Compute TE_C1
    TY_X[2] = HYX1X2 - HYX1 - HYY1X1X2 + HYY1X1
    # Compute TE
    TY_X[1] = HX1X2 - HX1 - HX1Y1X2 + HX1Y1
    # Compute MI
    TY_X[0] = HX2 + HY1 - HY1X2

    # Compute H(Xn+1|Xn(k),Yn-1(l))
    TY_X[4] = HYX1X2 - HYX1

    # Compute H(Xn+1|Xn(k+1),Yn-1(l))
    TY_X[5] = HXX1YX2 - HXX1Y

    return TY_X

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