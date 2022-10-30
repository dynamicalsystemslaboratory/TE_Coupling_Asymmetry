# this function computes the entropy of N time series of length T

import numpy as np
    
def compute_entropy(ts): 
    T = ts.shape[0] # no. of rows
    l,s,dcount = np.unique(ts,axis=0,return_inverse=True,return_counts=True)
    l = l*1 # converts to 0,1
    p_xs = dcount/T
    # p_xs = np.histogram(s,bins=np.arange(0,len(l)))/T
    # p_xs = p_xs + (p_xs == 0)
    H = - sum(np.multiply(p_xs,np.log2(p_xs)))
    
    return H
