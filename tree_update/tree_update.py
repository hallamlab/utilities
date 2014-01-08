#!/usr/bin/python

# load some packages
from __future__ import division
try:
    import re
    import sys
    import os
    import types
    import argparse
    
    from libs.LCAComputation import *
    from libs.MeganTree import *
except:
    import_exception = """ Could not load some modules """
    print import_exception 
    sys.exit(3)
    
what_i_do = "I'm a script that does something"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('--ncbi_node', dest='ncbi_node', type=str, nargs=1,
                required=False, help='a selection of one or many input .fasta files (required)', default=None)
parser.add_argument('--ncbi_name', dest='ncbi_name', type=str, nargs=1,
               required=False, help='the name of the output frequeny matrix (required)', default=None)
parser.add_argument('--megan_lvl_map', dest='megan_lvl_map', type=str, nargs=1,
              required=False, help='the name of the output frequeny matrix (required)', default=None)
parser.add_argument('--megan_tree', dest='megan_tree', type=str, nargs=1,
              required=False, help='the name of the output frequeny matrix (required)', default=None)
parser.add_argument('--megan_map', dest='megan_map', type=str, nargs=1,
                required=False, help='the name of the output frequeny matrix (required)', default=None)

def process_ncbi(ncbi_node_file, ncbi_name_file, megan_lvl_map):
    
    handle = open(ncbi_name_file, 'r')
    lines = handle.readlines()
    handle.close()
    ncbi_names = {}
    for l in lines:
        fields = l.split("|")
        taxa_id = fields[0].strip('\t')
        name = fields[1].strip('\t')
        ncbi_names[taxa_id] = name
    
    # print out mapping file
    handle = open(megan_lvl_map, 'r')
    lines = handle.readlines()
    handle.close()
    
    megan_level_map = {}
    for l in lines:
        fields = [ x.strip() for x in l.rstrip('\n').split('\t')]
        megan_level_map[fields[1]] = fields[0]
        
    handle = open(ncbi_node_file, 'r')
    lines = handle.readlines()
    handle.close()
    
    # writeout format for LCAComputation & MeganTree
    outfile = open("temptree", 'w')
    megan_map_out = open("ncbi.map", 'w')
    for l in lines:
        fields = l.split("|")
        tax_id = fields[0].strip('\t')
        parent_id = fields[1].strip('\t')
        rank = fields[2].strip("\t")
        # kishori's pipeline format
        line = ncbi_names[tax_id] + "\t" + tax_id + "\t" + parent_id + "\n"
        if rank in megan_level_map:
            rank = megan_level_map[rank]
        else:
            rank = '0'
        megan_map_line = tax_id + "\t" + ncbi_names[tax_id] + "\t" + "-1" + "\t" + rank + "\n"
        outfile.write(line)
        megan_map_out.write(megan_map_line)
    outfile.close()
    megan_map_out.close()

# process megan_map functionality
def process_megan_tree(megan_tree_file, megan_map_file, nodes_to_add = None):
    print "in process_megan_tree"
    meganTree = None
    lca = LCAComputation(megan_tree_file)
    # 1. create tree structure from megan newick format
    
    # 2. add nodes to tree if possible
    
    # 3. write tree back out into newick format for megan
    
    
# the main function
def main():
    # parse arguments
    args = vars(parser.parse_args())
    
    if args['megan_tree'] != None and args['megan_map'] != None:
        process_megan_tree(args['megan_tree'][0], args['megan_map'])
    exit()
    # read NCBI node file (nodes.dmp)
    process_ncbi(args['ncbi_node'][0],args['ncbi_name'][0],args['megan_lvl_map'][0])
    
    # experiment with Kishori's megan tree
    meganTree = None
    lca = LCAComputation('temptree')
    # lca = LCAComputation('data/pipeline/ncbi_taxonomy_tree_small.txt')
    meganTree = MeganTree(lca)
    meganTree.build_tree() # insert all taxa into tree
    meganTree.printNewickTree('1') # create newick tree
    outtree = open('ncbi.tre', 'w')
    outtree.write(meganTree.output)
    outtree.close()
    
    exit()


# call the main function
if __name__ == "__main__":
   sys.exit(main())
