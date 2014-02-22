#!python


"""
process_id2plate3a.py

Created by Niels Hanson on 2013-08-16.
Copyright (c) 2013 Steven J. Hallam Laboratory. All rights reserved.
"""
from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright 2013"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Niels W Hanson"
__status__ = "Release"


try:
     import os
     import re
     import argparse
     import glob
     from os import makedirs, sys, remove
except:
    print """ Could not load some modules """
    print """ """
    sys.exit(3)

what_i_do = "Strip fasta file names"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_maps', type=str, nargs='+',
                required=True, help='a selection of one or many fosmid input .fasta files (required)', default=None)             
parser.add_argument('-o', dest='output_directory', type=str, nargs='?',
                required=True, help='the name a target output directory (required)', default=None)

def main(argv):
    args = vars(parser.parse_args())
    
    for f in args["input_maps"]:
        handle = open(f, "r")
        lines = handle.readlines()
        handle.close()
        out_handle = open(os.path.abspath(args["output_directory"]) + "/" + os.path.basename(f) + ".map.txt", "w")
        for l in lines:
            fields = l.split(" ")
            forward_backward = ""
            jgi_fb = fields[0].split(".")[1]
            if jgi_fb.find("b"):
                forward_backward = "CPF"
            else:
                forward_backward = "CPR"
            line =  fields[0] + "\t" + fields[0][0:4] + "-" + fields[7] + "." + forward_backward + "_" + fields[2]
            print line
            out_handle.write(line + "\n")
        out_handle.close()
    exit()


# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
