#!/usr/bin/env python

'''
Created on Jun 20, 2011

@author: jonathanfriedman
'''

from lib.SurveyMatrix import Survey_matrix as SM

    
def kwargs_callback(option, opt, value, parser,**kwargs):
    d = kwargs['d']
    d[option.dest] = value
    return d
  

def RunSparCC(counts_file, algo='SparCC', **kwargs):
    ## read counts data
    print 'reading data'
    temp   = SM()
    counts = temp.from_file(counts_file)
    
    ## Calculate correlations between components using SparCC
    print 'computing correlations'
    cor, cov = counts.basis_corr(algo=algo , **kwargs)
    
    ## write out results
    print 'writing results'
    cor_file = kwargs.get('cor_file', 'cor_mat_' + algo + '.out')
    cor.writetxt(cor_file)
    print 'wrote ' + cor_file
    if cov is not None:
        cov_file = kwargs.get('cov_file', 'cov_mat_' + algo + '.out')
        cov.writetxt(cov_file)
        print 'wrote ' + cov_file
    
    print 'Done!'
        

if __name__ == '__main__':
    ## parse input arguments
    from optparse import OptionParser
    kwargs = {}
    usage  = ('Compute the correlation between components (e.g. OTUs).\n' 
              'By default uses the SparCC algorithm to account for compositional effects.\n' 
              'Correlation and covariance (when applies) matrices are written out as txt files. \n'
              'Counts file needs to be a tab delimited text file where columns are samples and rows are components (e.g. OTUS).\n'
              ' See example/fake_data.txt for an example file.\n' 
              '\n'
              'Usage:   python SparCC.py counts_file [options]\n'
              'Example: python SparCC.py example/fake_data.txt -i 20 --cor_file=Cor_mat.out')
    parser = OptionParser(usage)
    parser.add_option("-c", "--cor_file", dest="cor_file", type = 'str',
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs}, 
                      help="File to which correlation matrix will be written.")
    parser.add_option("-v", "--cov_file", dest="cov_file", type = 'str',
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs},
                      help="File to which covariance matrix will be written.")
    parser.add_option("-a", "--algo", dest="algo", default='SparCC',
                      help="Name of algorithm used to compute correlations (SparCC (default) | pearson | spearman | kendall)")
    parser.add_option("-i", "--iter", dest = 'iter', type ='int', default=20,
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs},
                      help="Number of inference iterations to average over (20 default).")
    parser.add_option("-x", "--xiter", dest = 'xiter', type ='int', default=10,
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs},
                      help="Number of exclusion iterations to remove strongly correlated pairs (10 default).")
    parser.add_option("-t", "--thershold", dest = 'th', type ='float', default=0.1,
                      action="callback", callback= kwargs_callback, callback_kwargs = {'d':kwargs},
                      help= "Correlation strength exclusion threshold (0.1 default).")

    (options, args) = parser.parse_args()
    counts_file     = args[0]
    algo            = options.algo
    
    ## run SparCC    
    RunSparCC(counts_file, algo = algo, **kwargs)

    

