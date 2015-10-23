
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def PTDetect(x, E):

    # Local Variables: a, b, E, d, xL, i, P, T, x
    # Function calls: PTDetect, length
    #% Peak detection in data x for a given threshold E
    #%
    #% [P,T] = PTDetect(x, E)
    #%
    #% Jacobson, ML. Auto-threshold peak detection in physiological signals,
    #% 2001.
    #%
    #% compiled: jsalvi@rockefeller.edu
    P = np.array([])
    T = np.array([])
    a = 1.
    b = 1.
    i = 0.
    d = 0.
    xL = length(x)
    while i != xL:
        i = i+1.
        if d == 0.:
            if x[int(a)-1] >= x[int(i)-1]+E:
                d = 2.
            elif x[int(i)-1] >= x[int(b)-1]+E:
                d = 1.
                
            
            if x[int(a)-1]<=x[int(i)-1]:
                a = i
            elif x[int(i)-1]<=x[int(b)-1]:
                b = i
                
            
        elif d == 1.:
            if x[int(a)-1]<=x[int(i)-1]:
                a = i
            elif x[int(a)-1] >= x[int(i)-1]+E:
                P = np.array(np.hstack((P, a)))
                b = i
                d = 2.
                
            
            
        elif d == 2.:
            if x[int(i)-1]<=x[int(b)-1]:
                b = i
            elif x[int(i)-1] >= x[int(b)-1]+E:
                T = np.array(np.hstack((T, b)))
                a = i
                d = 1.
                
            
            
        
        
    return [P, T]