'''
Created on Jul 22, 2010

@author: jonathanfriedman
'''
import scipy
import pylab
import scipy.cluster.hierarchy as sch
import figure_tools as FT
from Bio.Cluster import distancematrix
from numpy import arange



def clust_data(data, row_metric='euclidean', col_metric='euclidean', file = None, **kwargs):
    '''
    Take data matrix and do hierarchical clustering of rows and cols.
    Plot sorted heatmap.
    '''
    n,m = data.shape
    
    ## parse input args
    if 'frame' not in kwargs: frame = True
    sort_rows  = kwargs.get('sort_rows', True)
    sort_cols  = kwargs.get('sort_cols', True)
    
    
    ## set figure spacing
    if 'row_labels' in kwargs: 
        row_labels = kwargs['row_labels']
        if 'row_label_width' in kwargs: row_label_width = kwargs['row_label_width']
        else:                           row_label_width = max(map(lambda s: len(s),row_labels))*.05     
    else:   row_label_width = .01
    if 'col_labels' in kwargs: 
        col_labels = kwargs['col_labels']
        if 'col_label_width' in kwargs: col_label_width = kwargs['col_label_width']
        else:                           col_label_width = max(map(lambda s: len(s),col_labels))*.05
    else: col_label_width = .01 
    edge_margin = .05
    cbar_width  = .02
    dend_width  = 0.15
    if sort_rows: data_width  = 1 - 2*edge_margin - row_label_width - dend_width - cbar_width - .03
    else:         data_width  = 1 - 2*edge_margin - row_label_width - cbar_width - .03
    if sort_cols: data_height = 1 - 2*edge_margin - col_label_width - dend_width - .03 
    else:         data_height = 1 - 2*edge_margin - col_label_width - .03 
    
    fig = pylab.figure(figsize=(8,8))
    # Compute and plot row dendrogram.
    if sort_rows:
        D_row = sch.distance.pdist(data, metric = row_metric) # row distance matrix
        D_row = sch.distance.squareform(D_row)
        drow_left   = edge_margin
        drow_bottom = edge_margin 
        drow_width  = dend_width
        drow_height = data_height
        ax1  = fig.add_axes([drow_left,drow_bottom,drow_width,drow_height], frame_on=frame)
        Y    = sch.linkage(D_row, method='average')
        Z1   = sch.dendrogram(Y, orientation='right')
        idx1 = Z1['leaves']
        ax1.set_xticks([])
        ax1.set_yticks([])
    else:
        drow_left   = edge_margin
        drow_bottom = edge_margin 
        drow_width  = 0
        drow_height = data_height
        idx1        = kwargs.get('row_order', arange(n)) 

    # Compute and plot col dendrogram.
    if sort_cols:
        D_col = sch.distance.pdist(data.transpose(), metric = col_metric) # row distance matrix
        D_col = sch.distance.squareform(D_col)
        dcol_left   = drow_left + drow_width + row_label_width 
        dcol_bottom = edge_margin + data_height + col_label_width
        dcol_width  = data_width
        dcol_height = dend_width
        ax2  = fig.add_axes([dcol_left,dcol_bottom,dcol_width,dcol_height], frame_on=frame)
        Y    = sch.linkage(D_col, method='average')
        Z2   = sch.dendrogram(Y, orientation='top')
        idx2 = Z2['leaves']
        ax2.set_xticks([])
        ax2.set_yticks([])
    else:
        dcol_left   = drow_left + drow_width + row_label_width 
        dcol_bottom = edge_margin + data_height + col_label_width
        dcol_width  = data_width
        dcol_height = 0
        idx2        = kwargs.get('col_order', arange(m))  


    # Plot sorted data matrix matrix.
    mat_left   = dcol_left 
    mat_bottom = drow_bottom
    mat_width  = data_width
    mat_height = data_height
    axmatrix   = fig.add_axes([mat_left,mat_bottom,mat_width,mat_height])
    plot_log   = kwargs.get('plot_log', False)
    if plot_log: data = pylab.log10(data)
    data = data[idx1,:]
    data = data[:,idx2]
#    im = axmatrix.pcolormesh(data, aspect='auto', origin='lower')
    im = axmatrix.matshow(data, interpolation = 'nearest', aspect='auto', origin='lower')

    ## plot labels
    if 'row_labels' not in kwargs: 
        axmatrix.set_yticks([])
        row_labels_sorted = None
    else:
        row_labels_sorted = map(lambda i:row_labels[i] ,idx1)
        row_labels_sorted.reverse()
        axmatrix.set_yticks(arange(len(row_labels_sorted)) + 0.0)
        axmatrix.set_yticklabels(row_labels_sorted) 
#        FT.format_ticks(axmatrix,xaxis = False)
    if 'col_labels' not in kwargs: 
        axmatrix.set_xticks([])
        col_labels_sorted = None
    else:
        col_labels_sorted = map(lambda i:col_labels[i] ,idx2)
        pylab.xticks(arange(len(col_labels_sorted))+.0, rotation = 90)
        xtickNames = pylab.setp(axmatrix, xticklabels=col_labels_sorted)

        
    # Plot colorbar.
    cbar_left   = mat_left + mat_width + .02
    cbar_bottom = mat_bottom
    cbar_height = data_width 
    axcolor = fig.add_axes([cbar_left,cbar_bottom,cbar_width,cbar_height])
    pylab.colorbar(im, cax=axcolor)
    if file is not None: fig.savefig(file)
    return data, row_labels_sorted, col_labels_sorted


def heatmap_clust(D, labels = None, label_width = None, frame = True, file = None):
    
    if labels is None: label_width = .01
    elif label_width is None: label_width = max(map(lambda s: len(s),labels))*.012
    edge_margin = .05
    cbar_width  = .02
    dend_width  = 0.15
    heat_width  = 1 - 2*edge_margin - label_width - dend_width - cbar_width - .03 
    
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
    drow_left   = edge_margin
    drow_bottom = edge_margin 
    drow_width  = dend_width
    drow_height = heat_width    
    ax1 = fig.add_axes([drow_left,drow_bottom,drow_width,drow_height], frame_on=frame)
    Y = sch.linkage(D, method='average')
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # Compute and plot second dendrogram.
    dcol_left   = drow_left   + drow_width + label_width 
    dcol_bottom = edge_margin + heat_width + label_width
    dcol_width  = heat_width
    dcol_height = dend_width
    ax2 = fig.add_axes([dcol_left,dcol_bottom,dcol_width,dcol_height], frame_on=frame)
    Y = sch.linkage(D, method='average')
    Z2 = sch.dendrogram(Y, orientation='top')
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # Plot distance matrix.
    mat_left   = dcol_left 
    mat_bottom = drow_bottom
    mat_width  = heat_width
    mat_height = heat_width
    axmatrix = fig.add_axes([mat_left,mat_bottom,mat_width,mat_height])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
#    im = axmatrix.pcolormesh(D, aspect='auto', origin='lower')
    im = axmatrix.matshow(D, interpolation = 'nearest', aspect='auto', origin='lower')
    
    
    ## plot labels
    if labels is None:
        axmatrix.set_xticks([])
        axmatrix.set_yticks([])
    else:
        row_labels_sorted = map(lambda i:labels[i] ,idx1)
        row_labels_sorted.reverse()
        axmatrix.set_yticks(arange(len(row_labels_sorted)) + 0.0)
        axmatrix.set_yticklabels(row_labels_sorted) 
#        FT.format_ticks(axmatrix,xaxis = False)

        col_labels_sorted = map(lambda i:labels[i] ,idx2)
        pylab.xticks(arange(len(col_labels_sorted))+.0, rotation = 90)
        xtickNames = pylab.setp(axmatrix, xticklabels=col_labels_sorted)
        

    # Plot colorbar.
    cbar_left   = mat_left + mat_width + .02
    cbar_bottom = mat_bottom
    cbar_height = heat_width 
    axcolor = fig.add_axes([cbar_left,cbar_bottom,cbar_width,cbar_height])
    pylab.colorbar(im, cax=axcolor)
#    pylab.show()
    if file is not None: fig.savefig(file)
    

if __name__ == '__main__':
    n = 100
    m = 110
    x = scipy.rand(n,m)
    row_inds = [1,3,5,7]
    col_inds = [0,2,4,6]
    x[1,:] *= 10
    x[row_inds,:] *= 10
    x[:,col_inds] *= 10
    metric='euclidean'
    
    row_labels = map(lambda i: str(i), range(n))
    col_labels = map(lambda i: str(i), range(m))
#    clust_data(x, metric)
    clust_data(x, metric, row_labels = row_labels, col_labels = col_labels)
    
    
    