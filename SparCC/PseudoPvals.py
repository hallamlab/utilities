#!/usr/bin/env python

'''
Created on May 31, 2011

@author: jonathanfriedman
'''
import cPickle as pickle
import numpy as np

from copy import deepcopy
from itertools import combinations
from lib.MatrixDictionary import MatrixDictionary as MD

def kwargs_callback(option, opt, value, parser,**kwargs):
    d = kwargs['d']
    d[option.dest] = value
    return d

def Run(cor_file, base_file, n, **kwargs):   
   ## load  'real' correlation
   temp = MD()
   cor  = temp.from_file(cor_file)
   ## load shuffled correlations  
   cors_sim  = []
   for i in xrange(n):
       file =  base_file + '_' + str(i) + '.txt'
       temp = MD()
       c    = temp.from_file(file)
       cors_sim.append(c)
   ## determine test type
   test_type = kwargs.get('type','two_sided')
   if test_type not in ('one_sided','two_sided'): raise ValueError('Unkown test type %s' %test_type)
   ## compute p-vals
   p_vals    = deepcopy(cor)
   for o1,o2 in combinations(p_vals.row_labels(), 2):
       vals   = np.array([c[o1][o2] for c in cors_sim])
       c_real = cor[o1][o2]
       if test_type == 'two_sided':
           n_sig = len(np.where(np.abs(vals)>np.abs(c_real))[0]) 
       elif test_type == 'one_sided':
           if c_real >= 0:  n_sig = len(np.where(vals>c_real)[0])
           else:            n_sig = len(np.where(vals<c_real)[0])
       p = 1.*n_sig/n
       p_vals[o1][o2] = p
       p_vals[o2][o1] = p
   ## write p-vals    
   out_file = kwargs.get('out_file', cor_file +'.pvals')
   p_vals.writetxt(out_file)    
   
   
if __name__ == '__main__':
    ## parse input arguments
    from optparse import OptionParser
    kwargs = {}
    usage  = ('Compute pseudo p-vals from a set correlations obtained from shuffled data.\n' 
              'Pseudo p-vals are the percentage of times a correlation at least as extreme as the "real" one was observed in simulated datasets. \n'
              'p-values can be either two-sided (considering only the correlation magnitude) or one-sided (accounting for the sign of correlations).\n'
              'Files containing the simulated correlations should be named [prefix]_[num].txt. The naming prefix is the second input argument, and the total number of simulated sets is the third.\n'
              '\n'
              'Usage:   python PseudoPvals.py real_cor_file sim_cor_file_prefix num_simulations [options]\n'
              'Example: python PseudoPvals.py example/basis_corr/cor_mat_sparcc.out example/pvals/sim_cor 5 -o pvals.txt -t one_sided')
    parser = OptionParser(usage)
    parser.add_option("-t", "--type", dest="type", type = 'str',
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs}, 
                      help="Type of p-values to computed.  oned-sided | two-sided (default).")
    parser.add_option("-o", "--out_file", dest="out_file", type = 'str',
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs}, 
                      help="Name of file to which p-values will be written.")
    (options, args) = parser.parse_args()
    real_cor_file   = args[0]
    sim_cor_file    = args[1]
    n               = int(args[2])
     
    ## write sample distance
    Run(real_cor_file, sim_cor_file, n, **kwargs)   
   
    