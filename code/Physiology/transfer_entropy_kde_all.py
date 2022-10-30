## Computes MI,TE,TE_C1,TE_C2 using kernel density
import numpy as np
from compute_entropy import *
from sklearn.neighbors import KernelDensity
import time

def transfer_entropy_kde_all(Xall,Yall,r,ktype):  
    # r: bandwidth
    # ktype: type of kernel (e.g. 'gaussian','tophat')
    # xnbins, ynbins: # of equally spaced points where we want to determine PDF to compute TE

    TY_X = np.zeros(4)  # MI, TE, TE*, TE**
 
    # This func calculates TE_Y_to_X = I( X_k : Y_k-1 | X_k-1 )
    # --> mutual information between X and past Y given past X is known
    # by calculating the joint probability of the three quantities through binning

    # X_(k-1) = X, X_(k) = X1, X_(k+1) = X2
    # Y_(k-1) = Y, Y_(k) = Y1, Y_(k+1) = Y2
    N = len(Xall)
    X = Xall[0:N-2].reshape((N-2,1))
    Y = Yall[0:N-2].reshape((N-2,1))
    X1 = Xall[1:N-1].reshape((N-2,1))
    Y1 = Yall[1:N-1].reshape((N-2,1))
    X2 = Xall[2:N].reshape((N-2,1))

    kde = KernelDensity(kernel=ktype, bandwidth=r)

    # Compute TE_C2
    DATA = np.concatenate((X,Y,X1,Y1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gXYX1Y1X2 = np.exp(log_pdf)

    DATA = np.concatenate((X,Y,X1,Y1),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gXYX1Y1 = np.exp(log_pdf)

    DATA = np.concatenate((X,Y,X1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gXYX1X2 = np.exp(log_pdf)
    
    DATA = np.concatenate((X,Y,X1),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gXYX1 = np.exp(log_pdf)
    
    TY_X[3] = (1/(N-2))*np.sum(np.log2((gXYX1Y1X2*gXYX1)/(gXYX1Y1*gXYX1X2)))

    # Compute TE_C1
    DATA = np.concatenate((Y,X1,Y1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gYX1Y1X2 = np.exp(log_pdf)

    DATA = np.concatenate((Y,X1,Y1),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gYX1Y1 = np.exp(log_pdf)

    DATA = np.concatenate((Y,X1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gYX1X2 = np.exp(log_pdf)
    
    DATA = np.concatenate((Y,X1),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gYX1 = np.exp(log_pdf)
    
    TY_X[2] = (1/(N-2))*np.sum(np.log2((gYX1Y1X2*gYX1)/(gYX1Y1*gYX1X2)))

    # Compute TE
    DATA = np.concatenate((X1,Y1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gX1Y1X2 = np.exp(log_pdf)

    DATA = np.concatenate((X1,Y1),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gX1Y1 = np.exp(log_pdf)

    DATA = np.concatenate((X1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gX1X2 = np.exp(log_pdf)
    
    DATA = X1
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gX1 = np.exp(log_pdf)
    
    TY_X[1] = (1/(N-2))*np.sum(np.log2((gX1Y1X2*gX1)/(gX1Y1*gX1X2)))

    # Compute MI
    DATA = np.concatenate((Y1,X2),axis=1)
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gY1X2 = np.exp(log_pdf)

    DATA = Y1
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gY1 = np.exp(log_pdf)

    DATA = X2
    kde.fit(DATA)
    log_pdf = kde.score_samples(DATA)
    gX2 = np.exp(log_pdf)
    
    TY_X[0] = (1/(N-2))*np.sum(np.log2((gY1X2)/(gY1*gX2)))

    if (max(1.0*(TY_X<0))>0):
        print("Error! TY_X = ",TY_X)

    return TY_X