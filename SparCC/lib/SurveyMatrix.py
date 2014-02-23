'''
Created on Oct 6, 2010

@author: jonathanfriedman
'''

import cPickle as pickle
import numpy as np
from MatrixDictionary import MatrixDictionary
from copy import deepcopy

class Survey_matrix(MatrixDictionary):
    '''
    Class for storing Survey counts (e.g. 16S survey counts).
    Rows (outer dictionary) correspond to OTUs, cols (inner dict) correspond to samples.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        
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
        
        
    def from_file(self, file, trans = False):
        '''
        Create object from data file.
        File header is sample ids, each row starts with otu id.
        '''
        f       = open(file,'rU')
        header  = f.readline()
        cols    = header.strip().split('\t')[1:]
        rows    = []
        counts  = []
        for line in f:
            fields = line.strip().split('\t')
            rows.append(fields[0])
            try:
                counts.append( map(lambda x: float(x), fields[1:]) )
            except:
                print line
                print fields
        f.close()
        
        self.from_matrix( np.array(counts), row_labels = rows, col_labels = cols)
        
        if trans: mat =  self.transpose()
        else:     mat = self
        return mat
    
    
    def total_counts(self, rows = False):
        '''
        Return a dictionary of total counts per row (otu) or col (sample)
        '''        
        mat, row_labels, col_labels = self.to_matrix()
        
        if rows: 
            axis   = 1
            labels = row_labels
        else:   
            axis = 0
            labels = col_labels
        
        tot = mat.sum(axis = axis)
        
        ## make total dictionary
        total = {}
        for label, t in zip(labels,tot): total[label] = t
        return total
        
    
    
    def filter_otus(self, min_avg = 0, min_avg_frac = 0, min_present = 0):
        '''
        Remove otu that don't meet filtering criteria.
        Returns new instance.
        '''
        n          = len(self.col_labels()) # number of samples
        otus_2_del = set([])                # init otus to delete 
        if min_avg: # filter otus with < min_avg counts on average
            min_tot      = n*min_avg
            total_counts = self.total_counts(rows = True)
            t            = filter(lambda otu: total_counts[otu] < min_tot, total_counts.keys()) 
            otus_2_del.update(set(t))
        if min_avg_frac: # filter otus with < min_avg_frac fraction on average
            fracs        = self.to_fractions('normalize')
            min_tot      = n*min_avg_frac
            total_counts = fracs.total_counts(rows = True)
            t            = filter(lambda otu: total_counts[otu] < min_tot, total_counts.keys()) 
            otus_2_del.update(set(t))
        if min_present: # filter otus that are present in < min_present samples
            mat_d, samples = self.to_dict()
            t             = filter(lambda otu: len( np.where(mat_d[otu]>0)[0] ) < min_present, mat_d.keys()) 
        mat_filtered = self.remove_rows(otus_2_del)
        return mat_filtered
    
    
    def filter_samples(self, min_tot):
        '''
        Remove samples that don't meet filtering criteria.
        Returns new instance.
        '''
        total_counts = self.total_counts(rows = False)
        otus_2_del   = filter(lambda otu: total_counts[otu] < min_tot, total_counts.keys()) 
        mat_filtered = self.remove_cols(otus_2_del)
        return mat_filtered
    
    def combine(self, other, default = 0):
        '''
        Combine 2 HMPmats to a single one.
        OTUs that are missing from one HMPmat will be added to all samples with value given by 'default' input param. 
        '''
        self_otus  = set(self.row_labels())
        other_otus = set(other.row_labels())
        self_new   = self.add_rows(other_otus - self_otus,  vals = [dict.fromkeys(self.col_labels(),default)])
        other_new  = other.add_rows(self_otus  - other_otus, vals = [dict.fromkeys(other.col_labels(),default)])
        
        
        other_new_t   = other_new.transpose()
        combined_t    = (deepcopy(self_new)).transpose()
        for sample, vals in other_new_t.iteritems():
            combined_t[sample] = other_new_t[sample]
        combined = combined_t.transpose()
        return combined
    
    
    def normalize(self):
        '''
        Normalize col counts by total.
        Returns new instance.
        '''
        mat, row_labels, col_labels = self.to_matrix()
        total_counts                = self.total_counts(rows = False)
        
        for i,col in enumerate(col_labels):
            mat[:,i] /= total_counts[col]
            
        normed = self.remove_rows(self.row_labels())
        normed.from_matrix(mat, row_labels, col_labels)
        return normed

        
    
    def to_fractions(self, method = 'dirichlet', **kwargs):
        '''
        Convert counts to fraction, either by simple normalization or adding pseudo counts, or dirichlet sampling.
        If dirichlet sampling is used, for each sample (col) fit a dirichlet distribution and sample the fraction from it.
        The prior is a uniform dirichlet (a = ones(len(otus)) )
        Return a new instance.
        '''
        if method is 'normalize': fracs = self.normalize()
        elif method is 'pseudo': 
            p_counts = kwargs.get('p_counts',1) 
            fracs = (self+p_counts).normalize()
        else:
            from numpy.random.mtrand import dirichlet
            mat, row_labels, col_labels = self.to_matrix()
            for i in range(len(col_labels)): # for each sample
                N        = mat[:,i]     # counts of each otu in sample
                a        = N+1          # dirichlet parameters
                f        = dirichlet(a) # fractions are random sample from dirichlet
                mat[:,i] = f
            fracs = self.remove_rows(self.row_labels())
            fracs.from_matrix(mat, row_labels, col_labels)
        return fracs
                    
    
    def avg_abund(self, otus):
        '''
        Get the average abundance of given outs.
        Return dictionary keyed by out id with values = average abundance across all samples.
        '''
        d, samples = self.to_dict()
        avg        = {}
        for otu in otus: 
            if otu in d: avg[otu] = (d[otu]).mean() 
            else:        avg[otu] = 0
        return avg
        
    
    
    def add_internal_node_counts(self,tree):
        '''
        Add rows corresponding to all internal nodes in tree.
        Counts for an internal node are the sum of counts of all its kids.
        Returns new instance.
        '''
        d, col_labels  = self.to_dict() # convert self to dictionary
        internal_nodes = filter(lambda k: not tree[k].is_leaf() ,tree.keys()) # names of all internal nodes
        for node in internal_nodes:
            kids   = node.split('-')
            counts = np.zeros(len(col_labels)) # init counts of internal node
            for kid in kids: # sum up counts of kids of nodes, if present in HMP matrix. 
                if kid in d: counts += d[kid] 
            d[node] = counts
        new_mat = self.remove_rows(self.row_labels())
        new_mat.from_dict(d, col_labels)
        return new_mat

      
    def discriminating_OTUs(self,other):
        '''
        For each OTU compute the difference between the distribution in self samples and other samples
        '''
        from scipy.stats import mannwhitneyu, ks_2samp
        self_mat,  row_labels, col_labels = self.to_matrix()
        other_mat, row_labels, col_labels = other.to_matrix()
        stat      = {} # dictionary of value of statistic for each row
        p_val     = {} # dictionary of p-value of statistic for each row
        mean_self  = {} # dictionary of means 
        mean_other = {} # dictionary of means 
        for i,row in enumerate(row_labels):
            s, p = ks_2samp(self_mat[i,:],other_mat[i,:])
            stat[row]       = s
            p_val[row]      = p
            a               = np.mean(self_mat[i,:])
            b               = np.mean(other_mat[i,:])
            mean_self[row]  = a
            mean_other[row] = b
        return stat, p_val, mean_self, mean_other
    
    
    def entropy(self, method = 'ML'):
        '''
        Calculate the entropy of each column using the R entropy package.
        Return a dictionary of keyed by col label.
        Note that different methods work with normalized/unnormalized data!
        '''
        from R_utilities import entropy
        H = {}
        mat , row_labels, col_labels = self.to_matrix()
        for c,x in zip(col_labels,mat.transpose()): H[c] = entropy(x, method = method)[0]
        return H
    
    
    def diversity(self, index = 'hill_1', method = 'ML'):
        '''
        Return the effective number of species per sample of order q (See Jost, 2008?).
        q = 0 => Richness (in this case a better estimator, such as, Chao1, is required).
        q = 1 => exp(entropy)
        q = 2 => 1/Simpson's index.
        Output is a dict, keyed by sample id.
        Note that some methods work only with unnormalized data!
        '''
        from survey.diversity import diversity
        n_eff = {}
        mat , row_labels, col_labels = self.to_matrix()
        for c,x in zip(col_labels,mat.transpose()):
            y = (np.array([x])).transpose() 
            n_eff[c] = diversity(y, [index], methods = [method])[0]    
        return n_eff
    
    
    def plot_diversity(self, indices, **kwargs):
        '''
        Plot the diversity of various orders of the samples.
        '''
        import matplotlib.pyplot as plt
        samples       = kwargs.get('samples', self.col_labels())
        sort_samples  = kwargs.get('sort_samples', True)
        sort_index    = kwargs.get('sort_index', 0) # index of q by which to sort
        methods       = kwargs.get('methods', ['ML']*len(indices))
        y_log         = kwargs.get('y_log', True)
        x_log         = kwargs.get('x_log', False)
        fs            = kwargs.get('fs', 16)
        sample_labels = kwargs.get('sample_labels', False) # index of q by which to sort
        label_rot     = kwargs.get('label_rot', 90)
        show          = kwargs.get('show', False)
        file          = kwargs.get('file', False)
        plt.figure()
        
        ## get the diversity of all orders
        nf = []
        for index,method in zip(indices,methods):
            nf.append( self.diversity(index, method = method) )
        
        ## determine order of samples for plotting
        if sort_samples:  
            samples        = self.col_labels()
            samples_sorted = sorted(samples, key = lambda x: nf[sort_index][x], reverse = True)
        else: 
            samples_sorted = samples
        
        for i,index in enumerate(indices):
            x = []
            y = []     
            for j,s in enumerate(samples_sorted):
                yt = nf[i][s]
                if y_log: yt = np.log10(yt) 
                y.append(yt)
                xt = j+1
                if x_log:  xt = np.log10(xt)
                x.append(xt)      
            plt.plot(x, y, '--o', label = index)

        if y_log: ylab = 'log10(Diversity)'
        else:     ylab = 'Diversity'
        if x_log: xlab = 'log10(Rank)'
        else:     xlab = 'Rank'
        plt.xlabel(xlab, fontsize = fs)
        plt.ylabel(ylab, fontsize = fs)
        if sample_labels: plt.xticks(np.arange(1, len(samples)+1), samples_sorted, rotation= label_rot)
        else:             plt.xticks(np.arange(1, len(samples)+1), rotation = label_rot)
        plt.legend()
        if file: plt.savefig(file)
        if show: plt.show()
    
    
    def scatter_diversity(self,q1,q2, method1 = 'ML', method2 ='ML', **kwargs):
        '''
        Make a scatter plot of the diversity of order q1 vs. diversity of order q2.
        '''
        import matplotlib.pyplot as plt
        samples       = kwargs.get('samples', self.col_labels())
        fs            = kwargs.get('fs', 16)
        show          = kwargs.get('show', False)
        file          = kwargs.get('file', False)
        plt.figure()
        
        ## get the diversities
        nf1 = self.diversity(q1, method = method1)
        nf2 = self.diversity(q2, method = method2)
        x = []
        y = []
        for s,n1 in nf1.iteritems():
            x.append(n1)
            y.append(nf2[s])
        
        ## make scatter plot
        plt.plot(x,y, 'o')
        xlab = 'Diversity order %.1f' %q1
        ylab = 'Diversity order %.1f' %q2
        plt.xlabel(xlab, fontsize = fs)
        plt.ylabel(ylab, fontsize = fs)
        plt.grid()
        
        ## add text box with spearman correlation value
        show_cor  = kwargs.get('show_cor', True)
        if show_cor:
            import scipy.stats as stats 
            r, p  = stats.spearmanr(x,y)  
            s     = 'spearman r = %.1f \np_val = %1.0e' %(r,p) 
            str_x = .05
            str_y = .95
            ax    = plt.gca()
            plt.text(str_x, str_y, s, bbox=dict(boxstyle="round", fc="0.8"), transform = ax.transAxes, ha = 'left',va = 'top')
        
        if file: plt.savefig(file)
        if show: plt.show()
        
        
    def basis_corr(self, algo='SparCC', **kwargs):
        '''
        '''
        import basis_correlations as basecor 
        counts_t, otus, samples = self.to_matrix()
        counts = counts_t.transpose()
        cor_med, cov_med = basecor.main(counts, algo=algo, **kwargs)
        cor  = MatrixDictionary()
        
        cor.from_matrix(cor_med, otus, otus)
        if cov_med is None:
            cov = None
        else:
            cov  = MatrixDictionary()
            cov.from_matrix(cov_med, otus, otus)  
        return cor, cov
        
    
    def CSI_test(self):
        '''
        '''
        fracs         = self.to_fractions()
        f, rows, cols = fracs.to_compositions() 
        p, V_base     = f.CSI_test()
        return p
    
    
    def get_pvals(self, **kwargs):
        '''
        make a p-value matrix for all pairs of OTUs 
        '''
        import HMPStructures.OTUnetwork as ONet
        from itertools import combinations
        method = kwargs.get('method','csi')
        iter   = kwargs.get('iter',20)
        th     = kwargs.get('th',1e-1)
        net    = kwargs.get('net',None)
        
        pvals = deepcopy(self)
        otus  = self.row_labels(
                                )
        ## make network & get uncorrelated otus
        for cutoff in np.linspace(.01,.2,20):
            p_vec = np.zeros(iter)
            ONet.remove_edges(net)
            ONet.add_edges(net, cutoff = cutoff, statistic = 'matrix', method = 'abs_larger') # statistic can be p_vals, q_vals or matrix
            ## get uncorrelated otus
            degree    = net.degree()
            nodes     = degree.keys()
            nodes0    = filter(lambda node: degree[node] ==0, nodes)
            if len(nodes0) < 6: continue
            remove  = filter(lambda otu: otu not in set(nodes0), self.keys())
            abunds0 = self.remove_rows(remove)
            for i in range(iter): p_vec[i] = abunds0.CSI_test()
            p = p_vec.mean()
            break
        ## go over all edges and asses significance
        for pair in combinations(otus,2):
            p1_vec   = np.zeros(iter)
            p2_vec   = np.zeros(iter)
            p12_vec  = np.zeros(iter)
            o1,o2    = pair
            
            remove0  = filter(lambda x: x in abunds0 ,pair)
            temp     = abunds0.remove_rows(remove0)
#            abunds1  = temp.add_rows(keys = [o1], vals= [self[o1]])
#            abunds2  = temp.add_rows(keys = [o2], vals= [self[o2]])
            abunds12 = temp.add_rows(keys = [o1,o2], vals= [self[o1],self[o2]])
            for i in range(iter): 
#                p1_vec[i]  = abunds1.CSI_test()
#                p2_vec[i]  = abunds2.CSI_test()
                p12_vec[i] = abunds12.CSI_test()
#            p1  = p1_vec.mean()
#            p2  = p2_vec.mean()
            p12 = p12_vec.mean()
            
            pvals[o1][o2] = p12
            pvals[o2][o1] = p12
        return pvals
            
    
    
    def make_comp_network(self, algo = 'crude'):
        '''
        Plot network from corresponding compositional data
        '''
        f, rows, cols = self.to_compositions()
        f2 = f.replace_zeros()
        net = f2.make_network(algo = algo, row_labels = rows, col_labels = cols)
        return net
   
   
    def sample_clustering(self, metric = 'JSsqrt', iter = 100, **kwargs):
        '''
        Do hierarchical clustering of samples.
        Clustering confidence is estimated by converting to fractions using dirichlet sampling multiple times.
        Each time a new dendrogram is written in newick format.
        All dendrograms are than summarized using a consensus tree. 
        
        '''
        from phylogenetics.phylip.dist_2_newick import write_dist_mat, neighbor
        from phylogenetics.phylip.consensus import consense, write_bootstraps
        import subprocess
        
        remove_files      = kwargs.get('remove_files', True)
        
        path              = kwargs.get('path', '')
        dist_base         = kwargs.get('dist_base', path + 'dist_mat_' + metric + '_')
        tree_base         = kwargs.get('tree_base', path + 'tree_' + metric + '_')
        neighbor_out_file = kwargs.get('neighbor_out_file', path + 'neighbor.out' )
        
        trees_file    = kwargs.get('trees_file', path + 'trees_' + metric + '.txt')
        trees_f       = open(trees_file,'wb')
        
        frac_method = kwargs.get('frac_method', 'dirichlet')
        fracs       = self.to_fractions()
        comps, otus, samples = fracs.to_compositions()
        if 'label_d' in kwargs: labels = map(lambda s:kwargs['label_d'][s] ,samples)
        else: 
            label_d = {}
            labels   = []
            for j, s in enumerate(samples):
                labels.append(str(j))
                label_d[s] = str(j) 
        
        for i in xrange(iter):
            print i
            ## convert to compositions
            fracs = self.to_fractions(frac_method)
            comps, otus, samples = fracs.to_compositions()
            
            ## compute and write distance matrix
            dist_file = dist_base + str(i) + '.txt'
            dist   = comps.distance(comps, metric = metric)
#            write_dist_mat(dist, labels , dist_file)
            write_dist_mat(dist, labels , dist_file)

            
            ## compute and write tree
            tree_file = tree_base + str(i) + '.txt'
            neighbor(dist_file, neighbor_out_file, tree_file, path = path)
            
            ## add tree to trees file
            f    = open(tree_file,'r')
            tree = f.readlines()
            f.close()
            for line in tree: trees_f.write(line) 
            
            ## remove dist and tree files
            if remove_files:
                cmd = 'rm ' + dist_file 
                subprocess.call(cmd,shell=True)
                cmd = 'rm ' + tree_file 
                subprocess.call(cmd,shell=True)
                
              
        trees_f.close()
        
        ## compute and write consensus tree
        con_out_file   = kwargs.get('con_out_file',   path + 'consensus_' + metric + '.out')
        con_tree_file  = kwargs.get('con_tree_file',  path + 'consensus_tree_' + metric + '.newick')
        boot_tree_file = kwargs.get('boot_tree_file', path + 'consensus_tree_boot_' + metric + '.newick')
        consense(trees_file, con_out_file, con_tree_file, path = path)
        write_bootstraps(con_tree_file, boot_tree_file) 
           
        return label_d
    
    
    def plot_rank_abundance(self, **kwargs):
        '''
        Make a rank abundance plot for samples.
        '''
        import matplotlib.pyplot as plt
        y_log    = kwargs.get('y_log', True)
        x_log    = kwargs.get('x_log', True)
        samples  = kwargs.get('samples', self.col_labels())
        fs       = kwargs.get('fs', 16)
        show     = kwargs.get('show', False)
        file     = kwargs.get('file', False)
        plt.figure()
        d, otus = self.transpose().to_dict()
        for s in samples:
            vals = d[s]
            vals_sort = sorted(vals[vals > 0], reverse = True)
            if y_log: y = np.log10(vals_sort) 
            else:     y = vals_sort
            x = np.arange(1, len(y)+1)
            if x_log: x = np.log10(x)
            plt.plot(x, y, '-x')
        plt.xlabel('Rank', fontsize = fs)
        if y_log: ylab = 'log10(Abundance)'
        else:     ylab = 'Abundance'
        if x_log: xlab = 'log10(Rank)'
        else:     xlab = 'Rank'
        plt.xlabel(xlab, fontsize = fs)
        plt.ylabel(ylab, fontsize = fs)
        if file: plt.savefig(file)
        if show: plt.show()
        
    
    def plot_range_abunance(self, **kwargs):
        '''
        Make a plot of OTU range vs average fraction.
        '''
        import matplotlib.pyplot as plt
        range_log = kwargs.get('range_log', False)
        frac_log  = kwargs.get('frac_log', True)
        otus      = kwargs.get('otus', self.row_labels())
        fs        = kwargs.get('fs', 16)
        show      = kwargs.get('show', False)
        file      = kwargs.get('file', False)
        
        plt.figure()
        
        ## get average fraction of each otu
        n_samples = len(self.col_labels()) 
        otu_avg = (self/n_samples).total_counts(rows = True)
        
        ## get range of each otu (= effective number of samples)
        method    = kwargs.get('method', 'ML')
        order     = kwargs.get('order', 1)
        otu_range = self.transpose().diversity(method = method, q = order)
        
        x = []
        y = []     
        for o in otus:
            yt = otu_range[o]
            if range_log: yt = np.log10(yt) 
            y.append(yt)
            xt = otu_avg[o]
            if frac_log:  xt = np.log10(xt)
            x.append(xt)      
        plt.plot(x, y, 'bo')
        
        if range_log: ylab = 'log10(Effective # samples)'
        else:         ylab = 'Effective # samples'
        if frac_log:  xlab = 'log10(<fraction>)'
        else:         xlab = '<fraction>'
        plt.xlabel(xlab, fontsize = fs)
        plt.ylabel(ylab, fontsize = fs)
        plt.grid()
        
        ## add text box with spearman correlation value
        show_cor  = kwargs.get('show_cor', True)
        if show_cor:
            import scipy.stats as stats 
            r, p  = stats.spearmanr(x,y)  
            s     = 'spearman r = %.1f \np_val = %1.0e' %(r,p) 
            str_x = .05
            str_y = .95
            ax    = plt.gca()
            plt.text(str_x, str_y, s, bbox=dict(boxstyle="round", fc="0.8"), transform = ax.transAxes, ha = 'left',va = 'top')
        if file: plt.savefig(file)
        if show: plt.show()
                
        