'''
Created on May 2, 2011

@author: jonathanfriedman

Distance functions between two vectors that are of general utility.
'''
import numpy as np
from numpy import array, log, log2, sum, isnan, all

def KLsym(x,y): #Symmetric KL divergence
    if not all((x==0) == (y==0)): raise ValueError('KL divergence not defined when only one distribution is 0.')
    x[x==0] = 1 # set values where both distributions are 0 to the same (positive) value. This will not contribute to the final distance.
    y[y==0] = 1
    d = 0.5*np.sum( (x-y)*(log2(x) - log2(y)) )
    return d
    
    
def JS(x,y): #Jensen-shannon divergence
    import warnings
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    d1 = x*log2(2*x/(x+y))
    d2 = y*log2(2*y/(x+y))
    d1[isnan(d1)] = 0
    d2[isnan(d2)] = 0
    d = 0.5*sum(d1+d2)    
    return d


def JSsqrt(x,y): 
    d = JS(x,y)
    return d**0.5

    
def Morisita(x,y):
    d = 1-2*sum(x*y)/(sum(x**2)+sum(y**2))
    return d

def chao_jaccard(x,y):
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr  # needed for loading R packages
    import rpy2.robjects.numpy2ri               # allows R functions to accept numpy arrays
    rfos = importr('fossil')
    d    = 1-array(rfos.chao_jaccard(x,y))
    return d[0]

def aitchison(x,y, center = 'mean'):
    lx = log(x)
    ly = log(y)
    if center is 'mean':     m = np.mean
    elif center is 'median': m = median
    clr_x = lx - m(lx)
    clr_y = ly - m(ly)
    d = ( sum((clr_x-clr_y)**2) )**0.5
    return d


def pdist(x, metric):
    import scipy.cluster.hierarchy as sch
    if metric in set(['KLsym', 'JS', 'JSsqrt', 'Morisita','chao_jaccard','aitchison']):
        if metric == 'KLsym':           f = KLsym
        elif metric == 'JS':            f = JS
        elif metric == 'JSsqrt':        f = JSsqrt 
        elif metric == 'Morisita':      f = Morisita
        elif metric == 'chao_jaccard':  f = chao_jaccard
        elif metric == 'aitchison':     f = aitchison
        else: raise ValueError('Unkown metric %s' %metric)
        D = sch.distance.pdist(x, f)
    else: 
        D = sch.distance.pdist(x, metric = metric) # row distance matrix
    D = sch.distance.squareform(D)
    return D


def cdist(x, y, metric):
    import scipy.cluster.hierarchy as sch
    if metric in set(['KLsym', 'JS', 'JSsqrt', 'Morisita','chao_jaccard','aitchison']):
        if metric is 'KLsym': f = KLsym
        elif metric is 'JS':  f = JS
        elif metric is 'JSsqrt': f = JSsqrt 
        elif metric is 'Morisita': f = Morisita
        elif metric is 'chao_jaccard':  f = chao_jaccard
        elif metric is 'aitchison':  f = aitchison
        D = sch.distance.cdist(x,y, f)
    else: 
        D = sch.distance.cdist(x,y, metric = metric) # row distance matrix
    return D

    