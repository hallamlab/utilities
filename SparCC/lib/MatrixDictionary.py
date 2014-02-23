'''
Created on Aug 6, 2010

@author: jonathanfriedman
'''

import cPickle as pickle
import numpy as np
from copy import deepcopy
from numpy import ndarray, matrix

class MatrixDictionary(dict):
    '''
    A class for a matrix represent by a nested dictionary.
    Each item in the outer dictionary is a row in the matrix. 
    Each row is itself a dictionary with keys = column headers, values = matrix val.
    
    All keys in the inner dictionary should be the same, as they represent the column headers! 
    '''

    def __init__(self, square = False, *args, **kwargs):
        '''
        Constructor
        '''
        self.square = square
#        dict.__init__(self, *args, **kwargs)


    def __add__(self, other):
        '''
        Element-wise addition.
        When both are MatrixDictionaries, only work if self and other have the same col and row headers.
        '''
        mat1, rows1, cols1 = self.to_matrix()
        if type(other) is MatrixDictionary:
            mat2, rows2, cols2 = other.to_matrix()
            if set(rows1) != set(rows2) or set(cols1) != set(cols2):
                raise IOError('Added objects must have the same row and col labels')
            mat = mat1 + mat2
        else:
            mat = mat1 + other
        new = self.remove_rows(self.row_labels())
        new.from_matrix(mat, rows1, cols1)
        return new
    
    def __sub__(self, other):
        '''
        Element-wise subtraction of two MatrixDictionaries.
        When both are MatrixDictionaries, only work if self and other have the same col and row headers.
        '''
        mat1, rows1, cols1 = self.to_matrix()
        if type(other) is MatrixDictionary:
            mat2, rows2, cols2 = other.to_matrix()
            if set(rows1) != set(rows2) or set(cols1) != set(cols2):
                raise IOError('Added objects must have the same row and col labels')
            mat = mat1 - mat2
        else:
            mat = mat1 - other
        new = self.remove_rows(self.row_labels())
        new.from_matrix(mat, rows1, cols1)
        return new
    
    def __mul__(self, other):
        '''
        Element-wise multiplication of two MatrixDictionaries.
        When both are MatrixDictionaries, only work if self and other have the same col and row headers.
        '''
        mat1, rows1, cols1 = self.to_matrix()
        if type(other) is MatrixDictionary:
            mat2, rows2, cols2 = other.to_matrix()
            if set(rows1) != set(rows2) or set(cols1) != set(cols2):
                raise IOError('Added objects must have the same row and col labels')
            mat = mat1 * mat2
        else:
            mat = mat1 * other
        new = self.remove_rows(self.row_labels())
        new.from_matrix(mat, rows1, cols1)
        return new
    
    def __div__(self, other):
        '''
        Element-wise division of two MatrixDictionaries.
        When both are MatrixDictionaries, only work if self and other have the same col and row headers.
        '''
        mat1, rows1, cols1 = self.to_matrix()
        if type(other) is MatrixDictionary:
            mat2, rows2, cols2 = other.to_matrix()
            if set(rows1) != set(rows2) or set(cols1) != set(cols2):
                raise IOError('Added objects must have the same row and col labels')
            mat = mat1 / mat2
        else:
            mat = mat1 / other
        new = self.remove_rows(self.row_labels())
        new.from_matrix(mat, rows1, cols1)
        return new
    
            
    def save(self,file):
        ''' pickles into file'''
        f=open(file,'w')
        pickle.dump(self,f)
        f.close()
        
    def load(self,file):
        ''' unpickles from file'''
        f=open(file,'r')
        temp=pickle.load(f)
        f.close()
        return temp
    
    def shape(self): return (len(self.row_labels()), len(self.col_labels()) )
    
    def row_labels(self): return self.keys()
    
    def col_labels(self): return self.values()[0].keys()
           
           
    def add_rows(self, keys, vals = [{}]):
        ''' 
        Add rows to self, with values = vals. 
        Return new instance.
        !!! Currently works only for non square matrix !!!
        Inputs:
            keys = [list] keys of new elements
            vals = [list]
                    if len(vals) = 1: make all values be vals[0]
                    else:             len(vals) must be eqaul len(keys), and each key gets the corresponding val.
        
        Note: The 'fromkeys' dict method could be used, but I had bad results with it assigning the same value to several keys,
              i.e., changing the values of one key would change the values of other keys!
        '''
        ## verify that col keys in all added rows are identical
        added_cols = set(vals[0].keys()) 
        for v in vals[1:]:
            if added_cols != set(v.keys()): 
                raise IOError('All added row dict must have the same keys, as they represent columns!')
        ## verify that col keys match the existing col keys
        if len(self):
            col_keys = set(self.col_labels())
            if col_keys != added_cols: raise IOError('Added row dict must have the same keys as existing rows, as they represent columns!')
        
        n       = len(keys)
        mat_new = deepcopy(self)
        if len(vals) == 1:   
            for key in keys: mat_new[key] = deepcopy(vals[0])
        elif len(vals) == n: 
            for key,val in zip(keys,vals): mat_new[key] = val
        else:
            raise IOError('vals list length must either be 1, or equal to the length of the keys list.')
        return mat_new 
    
    
    def add_cols(self, keys, vals = [{}]):
        '''
        Add rows to self, with values = vals. 
        Return new instance.
        Same args as add_rows
        '''
        trans     = self.transpose()
        new_trans = trans.add_rows(keys,vals)
        return new_trans.transpose()
    
    def remove_rows(self,labels):
        '''
        Remove all row with given labels.
        Return new instance.
        '''
        mat_new = deepcopy(self)
        for l in labels:  del mat_new[l]
        if hasattr(self, 'sqaure') and self.square: # also remove cols if square
            mat_new = mat_new.remove_cols(labels)   
        return mat_new
    
    
    def remove_cols(self, labels):     
        '''
        Remove all cols with given labels.
        Return new instance.
        If square matrix, use remove rows.
        '''
        trans     = self.transpose()
        trans_new = trans.remove_rows(labels)
        mat_new   = trans_new.transpose()
        return mat_new 
          
    
    def from_file(self, file, trans = False):
        '''
        Create object from data file.
        File header is sample ids, each row starts with otu id.
        '''
        f       = open(file,'r')
        header  = f.readline()
        cols    = header.strip().split('\t')[1:]
        rows    = []
        counts  = []
        for line in f:
            fields = line.strip().split('\t')
            rows.append(fields[0])
            counts.append( map(lambda x: float(x), fields[1:]) )
        f.close()
        
        self.from_matrix( np.array(counts), row_labels = rows, col_labels = cols)
        
        if trans: mat =  self.transpose()
        else:     mat = self
        return mat
    
      
      
    def from_matrix(self, mat, row_labels = None, col_labels = None):
        '''
        Creates MatrixDictionary from matrix, + row & col labels
        Inputs:
                mat = [np array] data matrix
                row_labels = [list] keys of rows
                col_labels = [list] keys of cols
        '''
        n,m = mat.shape
        if row_labels is None: row_labels = ['row_' + str(i) for i in xrange(n)]
        if col_labels is None: col_labels = ['col' + str(i) for i in xrange(m)]
        
        ## verify that matrix dimensions match labels dims
        if len(row_labels) != mat.shape[0] or len(col_labels) != mat.shape[1]:
            raise IOError('Matrix dimensions must match labels dimensions')
     
        for i, row in enumerate(row_labels):
            self[row] = {}
            for j, col in enumerate(col_labels):
                self[row][col] = mat[i,j]
                
    
    def from_dict(self, d, col_labels):
        '''
        Creates MatrixDictionary from dictionary + col labels
        Inputs:
            d = [dict] Each entry corresponds to a row. 
                       Keys are row labels and values are np arrays of values.
           col_labels = [list] keys of cols                
        '''
        for row, counts in d.iteritems():
            self[row] = {}
            for j, col in enumerate(col_labels):
                self[row][col] = counts[j]

                    
    def from_nested_dict(self, d):
        '''
        Convert regular nested dictionary to NestedDictionary object.
        '''
        keys = d.keys()
        self.add_rows(keys)
        for key,val in d.iteritems():
            self[key] = val
        
    
    def to_matrix(self):
        '''
        Return the corresponding np array + lists of row and col keys
        '''
        row_labels = sorted(self.row_labels())
        col_labels = sorted(self.col_labels())
        if hasattr(self, 'sqaure') and self.square: col_labels = row_labels ## check if this is a square matrix, and make sure that makes same order of rows and cols
        matrix = np.zeros((len(row_labels), len(col_labels) ))
        for i,row_label in enumerate(row_labels):
            for j, col_label in enumerate(col_labels):
                try:
                    matrix[i,j] = self[row_label][col_label]
                except:
                    a =3
        return matrix, row_labels, col_labels
    
    
    def to_dict(self):
        '''
        Return the corresponding dictionary + list of col labels. 
        Each entry corresponds to a row. 
        Keys are row labels and values are np arrays of values.
        '''
        d = {}
        matrix, row_labels, col_labels = self.to_matrix()
        for i, label in enumerate(row_labels): d[label] = matrix[i,:]
        return d, col_labels
        
        
    def to_list(self):
        '''
        Return the corresponding nested list + lists of row and col keys
        '''
        row_labels = sorted(self.row_labels())
        col_labels = sorted(self.col_labels())
#        if self.square: col_labels = row_labels ## check if this is a square matrix, and make sure that makes same order of rows and cols
        l = []
        for i,row_label in enumerate(row_labels):
            l.append([])
            for j, col_label in enumerate(col_labels):
                try:
                    l[i].append(self[row_label][col_label])
                except:
                    a =3
        return l, row_labels, col_labels
    
    def to_compositions(self):
        '''
        Return corresponding compositional data object
        '''
        from Compositions import CompData
        matrix, row_labels, col_labels = self.to_matrix()
        f_mat = CompData(matrix.transpose(), dtype=float)
        return f_mat, row_labels, col_labels 
        
        
                    
    def to_PairMatrix(self):
        '''
        Return corresponding PairMatrix object (soon to be obsolete).
        To be used with square matrices only!!!
        '''
        from HMPStructures.OTUnetwork import PairMatrix
        mat, row_labels, col_labels = self.to_matrix()
        PM  = PairMatrix(mat,  row_labels)
        return PM
    
        
    def vals_by_keys(self, key_pairs):
        '''
        Return a list of values corresponding to key_pairs.
        Inputs:
            key_pairs = [list] each element = [row_key, col_key].
        Outputs:
            vals = [list] values for each pair in key_pairs, in corresponding order.
        '''
        vals = map(lambda pair: self[pair[0]][pair[1]], key_pairs)
        return vals
        
        
                
    def transpose(self):
        '''
        Change the order of nesting of dictionaries.
        Same as transposing the matrix of data corresponding to the nested dictionary.
        Only works if all nested dictionaries have the same keys!
        '''
        in_keys  = self.values()[0].keys()
        
        trans    = self.remove_rows(self.row_labels())
        trans    = trans.add_rows(in_keys)
        
        for key_out, d_in in self.iteritems():
            for key_in, val in d_in.iteritems(): trans[key_in][key_out] = val
        
        return trans  
      
        
    def log10(self, eps = 1e-10):
        '''
        Return the log10 version of self.
        '''
        mat, row_labels, col_labels = self.to_matrix()
        logged_mat = np.log10(mat + eps)
        logged     = self.remove_rows(self.row_labels())
        logged.from_matrix(logged_mat, row_labels, col_labels)
        return logged
    
    def log(self, base = np.e, eps = 1e-10):
        '''
        Return the log (with given base) version of self.
        '''
        mat, row_labels, col_labels = self.to_matrix()
        logged_mat = np.log(mat + eps)/np.log(base)
        logged     = self.remove_rows(self.row_labels())
        logged.from_matrix(logged_mat, row_labels, col_labels)
        return logged
    
    
    def powered(self, base = np.e):
        '''
        Return base^self.
        '''
        mat, row_labels, col_labels = self.to_matrix()
        powered_mat = base**mat
        powered     = self.remove_rows(self.row_labels())
        powered.from_matrix(powered_mat, row_labels, col_labels)
        return powered
        
    
    def min(self):
        mat, row_labels, col_labels = self.to_matrix()
        return mat.min()
    
    def max(self):
        mat, row_labels, col_labels = self.to_matrix()
        return mat.max()
                  
    def writetxt(self, file, header_flag = True):
        '''
        Write self to txt file.
        '''
        f = open(file, 'w')
        cols        = self.col_labels()
        rows        = self.row_labels()
        cols_sorted = sorted(cols)
        rows_sorted = sorted(rows)
        lines = []
        if header_flag: f.write('\t'.join(['OTU_id'] + cols_sorted) + '\n')
        for row in rows_sorted:
            row_d = self[row]
            vals  = []
            for col in cols_sorted: vals.append('%.3f' %row_d[col])
            if header_flag: line = '\t'.join([row] + vals)
            else:           line = '\t'.join(vals)
            lines.append(line)
        f.write('\n'.join(lines))
        f.close()
        
    
    def writesorted(self, lineages, base_file, path, metric = 'euclidean'):
        '''
        Do hierarchical clustering of rows and cols and write the sorted data to files.
        '''
        from heatmap_clust import clust_data, heatmap_clust
        matrix, row_labels, col_labels = self.to_matrix()
        mat_sorted, row_labels_sorted, col_labels_sorted = clust_data(matrix, metric, row_labels = row_labels,col_labels = col_labels, row_label_width = .25)
        row_labels_sorted.reverse()
        
        mat_file = path + base_file + '_abunds_sorted.txt'
        np.savetxt(mat_file, mat_sorted)
        
        
        otu_file = path + base_file + '_OTUs_sorted.txt'
        f = open(otu_file,'w')
        for otu in row_labels_sorted:
            f.write('\t'.join([otu, lineages[otu]]) + '\n')
        f.close()
        
        samples_file = path + base_file + '_samples_sorted.txt'
        f = open(samples_file,'w')
        for sample in col_labels_sorted:
            f.write(sample + '\n')
        f.close()
        
        
        
    
        
    def plot(self, log = False, file = None, show = False):
        '''
        Plot the matrix with sorted rows and cols
        '''
        import matplotlib.pyplot as plt
        matrix, row_labels, col_labels = self.to_matrix()
        if    log: mat = np.log10(matrix)
        else:      mat = matrix
        fig = plt.figure()
        plt.imshow(mat,interpolation='nearest')
        ax = plt.gca()
        plt.colorbar()
        
        ## set labels
        plt.setp(ax,xticks = np.arange(len(col_labels)))
        plt.setp(ax,yticks = np.arange(len(row_labels)))
        xtickNames = plt.setp(ax, xticklabels=col_labels)
        ytickNames = plt.setp(ax, yticklabels=row_labels)
 
        if file: plt.savefig(fig)
        if show: plt.show()        
        
            
    def plot_heatmap(self, row_metric  ='euclidean', col_metric  ='euclidean', file = None, **kwargs):
        '''
        Plot heatmap of self, sorted by heirarchical clustering with given distance metric
        '''
        from heatmap_clust import clust_data, heatmap_clust
        matrix, row_labels, col_labels = self.to_matrix()
        plot_row_labels = kwargs.get('plot_row_labels',False)
        plot_col_labels = kwargs.get('plot_col_labels',False)
        if plot_row_labels and plot_col_labels: 
            clust_data(matrix, row_metric, col_metric, file = file, row_labels = row_labels, col_labels = col_labels,  **kwargs)
        elif plot_row_labels and not plot_col_labels:
            clust_data(matrix, row_metric, col_metric, file = file, row_labels = row_labels,  **kwargs)
        elif not plot_row_labels and plot_col_labels:
            clust_data(matrix, row_metric, col_metric, file = file, col_labels = col_labels,  **kwargs)    
        else:
            clust_data(matrix, row_metric, col_metric, file = file, **kwargs)
        
    
    
    def dist_mat(self, metric  ='euclidean', transpose = True):
        '''
        Return a MD square distance matrix corresponding to distance between rows
        '''
        import distances
        if transpose: data    = self.transpose()
        else:         data    = self.transpose()
        mat, row_labels, col_labels = data.to_matrix()
        D_mat = distances.pdist(mat, metric)
        D     = MatrixDictionary()
        D.from_matrix(D_mat, row_labels, row_labels)
        return D
        
        
    def plot_dist_mat(self, metric  ='euclidean', file = None, transpose = True, show_labels = False, **kwargs):
        import distances
        from heatmap_clust import clust_data, heatmap_clust
        matrix, row_labels, col_labels = self.to_matrix()
        if transpose: 
            mat    = matrix.transpose()
            labels = map(lambda s: s.split('_')[-1], col_labels)
        else:         
            mat    = matrix
            labels = map(lambda s: s.split('_')[-1], row_labels)
        D = distances.pdist(mat, metric)
        if show_labels: heatmap_clust(D, file = file, labels =labels, **kwargs)
        else:           heatmap_clust(D, file = file, **kwargs)
    
        
    def fuzzy_clustering(self, k, r = 2, metric = 'euclidean', rows = True):
        '''
        Perform fuzzy c-means clustering on rows or cols.
        k = number of clusters.
        r = fuzziness exponent. Less fuzzy as r -> 1.
        '''
        import scipy.cluster.hierarchy as sch
        from R_utilities import c_means
        import distances
        
        if rows: matrix, row_labels, col_labels = self.to_matrix()
        else:    matrix, row_labels, col_labels = self.transpose().to_matrix()
        
        D = distances.pdist(matrix, metric)
        
        ## cluster each row and return a dict of cluster membership
        memb, memb_hard, stats = c_means(D,k,r, diss = True)
        membership             = MatrixDictionary()
        membership.from_matrix(memb, row_labels, range(1,k+1))
        membership_hard = {}
        for row, m in zip(row_labels, memb_hard): membership_hard[row] = m
        return membership, membership_hard, stats
    
    
    def correlation(self, type = 'pearson'):
        '''
        Calculate the correlation between all rows.
        Return matrix-dicts of tau and p-val.
        '''
        import scipy.stats as stats
        if type not in set(['pearson','kendall','spearman']): raise IOError('Specified correlation type is not supported.')
        mat, row_labels,col_labels = self.to_matrix()
        c = MatrixDictionary(square = True)
        p = MatrixDictionary(square = True)
        for i in range(len(row_labels)):
            ri    = row_labels[i]
            if ri not in c:
                c[ri] = {}
                p[ri] = {}
            for j in range(i, len(row_labels)):
                rj = row_labels[j]
                if rj not in c:
                    c[rj] = {}
                    p[rj] = {}
                if i == j: 
                    c[ri][ri] = 1
                    p[ri][ri] = 1
                    continue
                if type   == 'pearson':  c_temp,p_temp = stats.pearsonr(mat[i,:], mat[j,:])
                elif type == 'kendall':  c_temp,p_temp = stats.kendalltau(mat[i,:], mat[j,:])
                elif type == 'spearman': c_temp,p_temp = stats.spearmanr(mat[i,:], mat[j,:])
                c[ri][rj] = c_temp
                c[rj][ri] = c_temp
                p[ri][rj] = p_temp
                p[rj][ri] = p_temp
        return c,p
    
    
    def MI(self, method = 'ML'):
        '''
        Compute Mutual Information, Variation of Information and Coefficient of Constraint between all rows.
        Currently only works for binary values.
        Return Matrix dictionary object of MIs. 
        '''
        import information_theory as IT
        l, row_labels, col_labels = self.to_list()
        info   = IT.Information()
        MI_mat,VI_mat,CC_mat = info.MI(l, method = method, alphabet = [[0,1]], shared_alph = True)
        MI = MatrixDictionary(square = True)
        VI = MatrixDictionary(square = True)
        CC = MatrixDictionary(square = True)
        MI.from_matrix(MI_mat, row_labels, row_labels)
        VI.from_matrix(VI_mat, row_labels, row_labels)
        CC.from_matrix(CC_mat, row_labels, row_labels)
        return MI,VI,CC
    
    
    def to_binary(self, th = 0):
        '''
        Discretize matrix s.t. matrix[matrix > th] = 1, matrix[matrix <= th] = 0. 
        Return new instance.
        '''
        mat, row_labels, col_labels = self.to_matrix()
        mat_bin = np.zeros(mat.shape)
        mat_bin[mat > th]  = 1
        bin     = self.remove_rows(self.row_labels())
        bin.from_matrix(mat_bin, row_labels, col_labels)
        return bin
    
    
    def sim_data(self, n, method = 'permute', base_file = None, format = 'txt'):
        '''
        Create n simulated data sets using method, having the same properties as the original.
        ('same' properties mean different things for different methods)
        If base_file is given, save each simulated data set with a running number starting from 0.
        Return a list of the simulatd data sets. 
        '''
        from simulate_data import permute_w_replacement as pwr
        if method is 'permute': sim_fun = pwr
        sims = []
        for i in xrange(n):
            sim = sim_fun(self)
            if base_file: # save simulated data
                if format =='txt':
                    file = base_file + '_' + str(i) + '.txt'
                    sim.writetxt(file)
                elif format == 'pick':
                    file = base_file + '_' + str(i) + '.pick'
                    sim.save(file)
                else:
                    raise ValueError('Unkown format %s' %format)              
            sims.append(sim)
        return sims
                    

def test_c_means():
    mat = np.array([[0.1,0],[0,0],[0,0],[0,0],[0,0],[1,1],[1,1],[1,1],[1,1],[1,1]])
    mat = np.random.rand(10,2)
    print mat
    col_labels = ['c0', 'c1']
    row_labels = range(10)
    data = MatrixDictionary()
    data.from_matrix(mat, row_labels, col_labels)
    metric = 'euclidean'
    memb, membership_hard, stats = data.fuzzy_clustering(2, metric = metric)
    print memb.to_matrix()
    print membership_hard
    print stats

def test_corr():
    n = 5
    m = 10
    mat = np.random.rand(n,m)
    mat[1,:] = mat[0,:]
    row_labels = range(n)
    col_labels = range(m)
    data = MatrixDictionary()
    data.from_matrix(mat, row_labels, col_labels)
    t,p = data.correlation(type = 'kendalltau')
    t_mat, row,col = t.to_matrix()
    plt.matshow(t_mat)
    plt.show()

def test_MI():
    n = 3
    m = 50
    mat = np.random.randint(0,2,(n,m))
    mat[2,:] = mat[0,:]
    row_labels = range(n)
    col_labels = range(m)
    data = MatrixDictionary()
    data.from_matrix(mat, row_labels, col_labels)
    print data.MI()
    
    
    
if __name__ == '__main__':
    
#    test_sim_data()
     
    row_labels = ['r0','r1','r2']
    col_labels = ['c0','c1']
    mat  = np.array([ [00,2],[1,11],[20,21] ])
    d = MatrixDictionary()
    d.from_matrix(mat, row_labels, col_labels)
    print d
    print d.sim_data(5)
#    print d.powered(base = 2) 
#    print d/mat
#    mat, rows, cols =  d.to_matrix()
#    print d.vals_by_keys([['r0','c0'],['r2','c1'] ])
#    bin = d.to_binary(th = 0)
     
#    d = MatrixDictionary()    
#    keys = ['sam1', 'sam2', 'sam3']
#    val  = {'par1': 0, 'par2' :0}
#    val2 = {'par1': 1, 'par2' :1}
#    d.add_rows(keys, [val,val,val2])
#    d['sam1']['par1'] = 11
#    d['sam1']['par2'] = 12
#    d['sam3']['par1'] = 31
#    print d.to_list()[0], '\n'
#    print d.transpose().to_list()[0] 
    
#    
#    path = '/Users/jonathanfriedman/Documents/Alm/HMP_HGT/blastn/'
#    f = open(path + 'hmp_pair_shared.pick')
#    shared = pickle.load(f)
#    s = MatrixDictionary() 
#    s.from_dict(shared)
#    s.writetxt(path + 'hmp_pair_shared.txt')
#    f.close()      
             
        
        
        
        