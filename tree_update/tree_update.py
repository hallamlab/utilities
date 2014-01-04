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
parser.add_argument('-m', dest='metapathways_ncbi_tree_file', type=str, nargs='+',
                required=True, help='a selection of one or many input .fasta files (required)', default=None)
# parser.add_argument('-o', dest='output_file', type=str, nargs='?',
#               required=True, help='the name of the output frequeny matrix (required)', default=None)


# the main function
def main():
    # parse arguments
    args = vars(parser.parse_args())
    
    # experiment with Kishori's megan tree
    meganTree = None
    lca = LCAComputation(args['metapathways_ncbi_tree_file'][0])
    meganTree = MeganTree(lca)
    
    exit()


# call the main function
if __name__ == "__main__":
   sys.exit(main())
