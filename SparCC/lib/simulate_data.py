'''
Created on Oct 15, 2010

@author: jonathanfriedman

Methods for creating fake HMP like data. 
'''
import numpy as np
import cPickle as pickle
from SurveyMatrix import Survey_matrix
#from estimate_dist import estimate
#from numpy.random.mtrand import dirichlet
#from numpy.random import multinomial, binomial


def independent_betabinom(original):
    counts , otus, samples = original.to_matrix()
    k = len(otus)
    n = len(samples)
    N = original.total_counts(rows = False).values()
    N_avg = np.mean(N)
    ## estimate beta-dist params for each OTU
    a = []
    tot = counts.sum(axis = 0)
    for i in range(k):
        c = np.c_[counts[i,:], tot- counts[i,:]]
        a_opt, ll = estimate.dirmulti(c, symmetric = False)
        a.append(a_opt)
    ## simulate each OTU
    simulated = np.zeros((k,n))
    for i in range(k):
        probs = dirichlet(a[i],n)
        for j,p in enumerate(probs):
            simulated[i,j] = multinomial(N_avg,p)[0]
    ## convert to Survey_matrix
    sim = Survey_matrix()
    sim.from_matrix(simulated, otus, samples)
    return sim
        

def independent_kde(site, original):
    from numpy.random import rand        
    counts , otus, samples = original.to_matrix()
    k = len(otus)
    n = len(samples)
#    N = original.total_counts(rows = False).values()
#    N_avg = np.mean(N)
    N = [1000] * n
    ## load kde of each OTU
    path  = '../data/otu_dist/'
    file  = site + '_filtered_kde.pick'
    f     = open(path + file,'r')
    kdes, cdfs  = pickle.load(f)
    f.close()
    ## simulate each OTU
    x = np.linspace(0,1,100)
    simulated = np.zeros((k,n))
    for i,otu in enumerate(otus):
        cdf = cdfs[otu]
        for j in range(n):
            p = x[cdf>rand()][0]  
            simulated[i,j] = binomial(N[j],p)
    ## convert to Survey_matrix
    sim = Survey_matrix()
    sim.from_matrix(simulated, otus, samples)
    return sim
        

def permute_w_replacement(original):
    '''
    Create simulated dataset were the count of each otu in each sample are randomly sampled from the all the counts of that otu in all samples.
    '''
    from numpy.random import randint 
    counts , otus, samples = original.to_matrix()
    n,k = counts.shape
    new = np.zeros( (n,k) )
    for row in range(n):
        new[row] = counts[row, randint(0,k,(1,k)) ]
    sim = Survey_matrix()
    sim.from_matrix(new, otus, samples)
    return sim
    
    


def dirichlet_multi(a , n, N, **kwargs):
    '''
    Each sample is a sample from a multinomial with probabilities drawn from a dirichlet.
    Inputs:
        a = [array] alpha parameters of dirichlet. length determins the number of otus.
        n = [int] number of samples.
        N = [int/array] number of reads to take for each sample.
            If int take same number of reads from each sample.
            If list, length must be = n.
    Optional inputs:
        otu_labels    = [list] row labels for Survey_matrix. Default to ['otu_1','otu_2',...,'otu_k'].
        sample_labels = [list] col labels for Survey_matrix. Default to ['sample_1','sample_2',...,'sample_n'].
    Return Survey_matrix of simulated data.
    '''    
    k     = len(a)              # number of otus
    if 'otu_labels' not in kwargs: 
        otu_labels   = map(lambda i: 'otu_' + str(i) ,range(1,1+k))
    else: otu_labels = kwargs['otu_labels']
    if 'sample_labels' not in kwargs: 
        sample_labels = map(lambda i: 'sample_' + str(i) ,range(1,1+n))
    else: sample_labels = kwargs['sample_labels']
    if type(N) == int: N = [N]*n
    
    probs = dirichlet(a,int(n)) # each row gives the probabilities of different otus in 1 sample.
    mat   = np.zeros((k,n))
    for i in range(n):
        mat[:,i] =  multinomial(N[i],probs[i,:])
    sim_data   = Survey_matrix()
    sim_data.from_matrix(mat ,row_labels = otu_labels, col_labels = sample_labels )
    return sim_data
    

def simulate(original, sim_fun, num_sim = 1e3, a = None, symmetric = True):
    '''
    Create num_sim simulated data sets using sim_fun, having the same properties as the original.
    ('same' properties mean different things for different sim_funs)
    '''
    from numpy import tile
    sims    = []
    counts , otus, samples = original.to_matrix()
    k = len(otus)
    n = len(samples)
    N = original.total_counts(rows = False).values()
    if a is None:
        a,ll  = estimate.dirmulti(counts.transpose(), symmetric = symmetric)
    if sim_fun ==  dirichlet_multi:
        for i in range(num_sim):
#            if not i%10: print i
            sims.append(dirichlet_multi(a , n, N,otu_labels=otus, sample_labels = samples))
    return sims
 
 
def test_dirichlet():
    import matplotlib.pyplot as plt
    N = 100
    n = 200
    k = 50
    a0 = .1
    a1 = .9*a0
    a = np.array([a1]+[(a0-a1)/(k-1)]*(k-1))
    abunds =  dirichlet_multi(a , n, N)
    abunds.plot_heatmap(metric  ='euclidean')
    plt.show()
    

def test_simulate():
    mat = Survey_matrix()
    x = np.array([[20,20],[20,40]])
    r_lab = ['r0','r1']
    c_lab = ['c0','c1']
    mat.from_matrix(x,r_lab,c_lab)
    sims = simulate(mat, dirichlet_multi, num_sim = 10)
    print sims


def test_permute():
    mat = Survey_matrix()
    x = np.array([[1,2,3,4,5],[6,7,8,9,10]])
    r_lab = ['r0','r1']
    c_lab = ['c0','c1','c2','c3','c4']
    mat.from_matrix(x,r_lab,c_lab)
    sims = permute_w_replacement(mat)
    print x
    print sims.to_matrix()
    

def test_independent_betabinom(): 
    n = 5
    m = 50
    mat = np.floor(100*np.random.rand(n,m))
    mat[1,:] = 100-mat[0,:]
    row_labels = range(n)
    col_labels = range(m)
    data = Survey_matrix()
    data.from_matrix(mat, row_labels, col_labels)
    c,p = data.correlation(type = 'pearsonr')
    print c.to_matrix()
    print p.to_matrix()
    sim = independent_betabinom(data)
    c_sim, p = sim.correlation(type = 'pearsonr')
    print c_sim.to_matrix()
    print p.to_matrix()
    sim_norm = sim.normalize(rows = False)
    c_sim_norm, p = sim_norm.correlation(type = 'pearsonr')
    print c_sim_norm.to_matrix()
    print p.to_matrix()
    
    
        
def main():
        n = 100
        method = 'permute'
#        sims = simulate(abunds, dirichlet_multi, num_sim = n)
        if method is 'betabinom': sim = independent_betabinom(abunds)
        elif method is 'kde':     sim = independent_kde(site, abunds)
        elif method is 'permute': sim = permute_w_replacement(abunds)
        sims = [sim]
#        ## normalize 
#        print '############ normalizing'
#        sims_norm = map(lambda s: s.normalize(rows = False), sims)
        ## save
        print '############ saving counts'
        if sparse_sim:  file = site +'_filtered_sim_lognorm_sparse_sim_' + method + '.pick'
        else:           file = site + '_filtered_sim_' + method + '.pick'
        f    = open(path + file,'w')
        pickle.dump(sims,f)
        f.close()
#        print '############ saving fractions'
#        file = site + '_filtered_sim_kde_norm.pick'
#        f    = open(path + file,'w')
#        pickle.dump(sims_norm,f)
#        f.close()
              
    

if __name__ == '__main__':
    main()
    
    
    