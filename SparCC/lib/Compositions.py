'''
Created on Nov 5, 2010

@author: jonathanfriedman
'''
import numpy as np
from copy import deepcopy
from numpy import array, zeros, tile, ones, log, log2, corrcoef, var, exp, cov, r_, diag, matrix, diag, sqrt, where, multiply, mean
from numpy.random.mtrand import multivariate_normal as mvn
from SurveyMatrix import Survey_matrix
from MatrixDictionary import MatrixDictionary

class CompData(matrix):
    '''
    Matrix of compositional data.
    Columns are components and rows are samples.
    When the data is closed, rows sum to a constant.
    '''

    def __init__(self, *args, **kwargs):
        '''
        Constructor
        '''
        basis = deepcopy(self)
        self.basis = basis
        self /= self.sum(axis = 1)
        
    def __add__(self, other):
        '''
        pertubation operation.
        If other is same size as self, each sample (row) in self will be perturbed by its counterpart in other.
        Elif, other a single sample, each sample in self is perturbed by the other sample.
        Else, error is raised.
        '''
        from numpy import multiply
        return CompData(multiply(self, other), dtype=float)
        
           
    def close(self):
        '''
        Normalize the row data.
        Return new instance. 
        '''
        closed = matrix(self, dtype=float)
        closed /= closed.sum(axis = 1)
        return closed


    def alr(self, comp = -1):
        '''
        Compute the additive log-ratio (alr) transformation of self with respect to the component given in comp.
        ''' 
        n,k = self.shape
        norm = tile(self[:,comp],(1,k-1))    
        y    = np.delete(self,comp,axis = 1)
        y /= norm
        return array(log(y))

    
    def clr(self, **kwargs):
        '''
        Do the central log-ratio (clr) transformation of self.
        kwarg 'centraility' is the metric of central tendency to divide by after taking the logarithm.
        '''
        centrality = kwargs.get('centraility' ,'mean')
        n,k  = self.shape
        temp = log(self)
        if   centrality is 'mean':   m = np.mean(temp,axis = 1)
        elif centrality is 'median': m = np.median(temp,axis = 1)
        z = temp - tile(m,(1,k))
        return z
    
    
    def variation_mat(self, shrink = False):
        '''
        Return the variation matrix of self.
        Element i,j is the variance of the log ratio of components i and j.
        Can also employ shrunken estimation of the variation matrix, following Rainer Opgen-Rhein and Korbinian Strimmer 2007a.
        '''
        self_a = array(self)
        n,k    = self.shape
        V      = zeros((k,k))
        y_mat  = []
        for i in range(k-1):
            for j in range(i+1,k):
                y     = array(log(self_a[:,i]/self_a[:,j]))
                y_mat.append(y)
                v = var(y, ddof = 1) # set ddof to divide by (n-1), rather than n, thus getting an unbiased estimator (rather than the ML one). 
                V[i,j] = v
                V[j,i] = v
        if not shrink: return matrix(V)
        else:
            from R_utilities import var_shrink
            data         = (array(y_mat)).transpose()
            V_shrink_vec = var_shrink(data)
            V_shrink     = zeros((k,k))
            m = 0
            for i in range(k-1):
                for j in range(i+1,k):
                   V_shrink[i,j] = V_shrink_vec[m]
                   V_shrink[j,i] = V_shrink[i,j]
                   m += 1 
            return matrix(V_shrink) 
    
    
    def replace_zeros(self, **kwargs):
        '''
        Replace the zeros by a small value by imputation.
        Return new object.
        
        Inputs:
            e = [float] fraction of minimal value to use as imputed value delta.
        '''
        type = kwargs.get('type','multiplicative')
        e   = kwargs.get('e',0.5)
        new = deepcopy(self)
        n,k = new.shape
        for i in range(n):
            vals         = new[i,:]                # fractions in current sample
            inds_z       = (vals ==0).nonzero()    # indices of zeros
            inds         = (vals > 0).nonzero()    # indices of no zeros
            delta        = e * np.min(vals[inds])  # imputed value for current sample
            vals[inds_z] = delta                   # replace zeros by imputed values
            if type is 'simple':
                vals /= np.sum(vals)
            elif type is 'multiplicative':
                vals[inds] *= (1-delta*len(inds_z))   
            new[i,:] = vals
        return new 
        
        
    def distance(self, other, metric = 'aitchison', **kwargs):
        '''
        Compute the distance matrix between self and other compositional matrix.
        Inputs:
            metric = [str] type of distance metric to use. 
        Output: 
            D = [array] distance matrix. d_ij = distance between self sample i (self[i,:]) and other sample j.
        '''
        from scipy.spatial.distance import cdist
        import distances
        n,k = self.shape
        m,k = other.shape
        D   = zeros((n,m))
        D   = distances.cdist(self,other,metric)
        return D
            
    
    def CSI_test(self):
        '''
        Do a test for complete subcompositional independence of all components (columns) of f. 
        '''
        from compositions.CSI_test import CSI_test
        p, V_base = CSI_test(self)
        return p, V_base
    
    
    def basis_corr(self, method = 'sparse', **kwargs):
        '''
        Estimate the correlation in the basis of self, assuming that the average correlation is small.
        '''
        from basis_correlations import basis_corr as corr
        V_base, Cor_base, Cov_base = corr(self, method, **kwargs)
        return V_base, Cor_base, Cov_base
    
    
    def to_SurveyMat(self, **kwargs):
        temp       = Survey_matrix()
        n, k       = self.shape
        row_labels = kwargs.get('row_labels', range(k))
        col_labels = kwargs.get('col_labels', range(n))
        temp.from_matrix(self.transpose(), row_labels, col_labels)
        return temp
    
    
    def make_network(self, algo = 'crude', **kwargs):
        '''
        '''
        import SurveyStructures.OTUnetwork as ONet
        n, k       = self.shape
        row_labels = kwargs.get('row_labels', range(k))
        col_labels = kwargs.get('col_labels', range(n))
        
        md = self.to_Surveymat(**kwargs)
        if   algo == 'crude':   mat,pval = (md.log10()).correlation(type = 'pearsonr')
        elif algo in set(['sparse', 'csi', 'clr']):   
            v_sparse, c_sparse, Cov_sparse = self.basis_corr(method = algo,**kwargs)
            mat   = MatrixDictionary()
            mat.from_matrix(c_sparse, row_labels, row_labels)
            pval = mat
        elif algo == 'basis':
            md     = Survey_matrix()
            n, k   = self.shape
            mat    = np.log(self.basis)
            md.from_matrix(mat.transpose(), row_labels, col_labels)
            mat,pval = md.correlation(type = 'pearsonr')
        elif algo in set(['pearsonr','kendalltau','spearmanr']):
            abunds = self.to_Surveymat()
            mat,pval = abunds.correlation(type = algo)

                
        ## make network
        if algo == 'MI': net = ONet.OTUnetwork(  abunds = abunds, lineages = lineages, matrix = mat.to_PairMatrix(), p_vals = {'direct':None}, algo = 'MI', lag = 0)
        else:  net = ONet.OTUnetwork(node_ids = mat.row_labels(), matrix = mat, p_vals = {'direct':pval}, algo = algo, lag = 0)
        return net
        
    
    def plot_network(self, algo = 'crude', file = None, **kwargs):
        import pylab
        import SurveyStructures.OTUnetwork as ONet
        cutoff = kwargs.get('cutoff',0.1)
        show   = kwargs.get('show',True)
        ## make network
        net = self.make_network(algo,**kwargs)
        ## set layout
        nodes = net.nodes()
        n     = len(nodes)
        t = np.linspace(0,2*np.pi,n+1)
        x_vec = np.cos(t[:-1])
        y_vec = np.sin(t[:-1])  
        pos = {}
        for node,x,y in zip(nodes,x_vec,y_vec): pos[node] = np.array([x,y])
        ## add edges
        ONet.remove_edges(net)
        ONet.add_edges(net, cutoff = cutoff, statistic = 'matrix', method = 'abs_larger') # statistic can be p_vals, q_vals or matrix
        ## plot
        ONet.plot_network(net, layout = 'circular', file = file, **kwargs)
        return net
                 