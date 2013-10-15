#!/usr/bin/python
"""
metapathways_last_to_megan.py

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
     import sys
     import argparse
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)


# usage: python metapathways_last_to_megan.py -i input_filename -o output_filename
what_i_do = "Parses blastout/blastout.parsed files and formats them as csv files for import into MEGAN"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_file', type=str, nargs='?',
                required=True, help='the input file to be parsed', default=None)                
parser.add_argument('-o', dest='output_csv', type=str, nargs='?',
                required=True, help='the output csv file to be created', default=None)


def main(argv):
    args = vars(parser.parse_args())
    
    input_filename = args['input_file']
    output_filename = args['output_csv']
    file_handle = open(input_filename, "r")

    lines = file_handle.readlines()

    file_handle.close()

    brackets_pattern = re.compile("\[(.*?)\]")

    output_handle = open(output_filename, "w")

    for l in lines:
    	fields = l.split("\t")
    	hits = brackets_pattern.search(fields[8])
    	if hits:
    	   read = fields[0]
    	   last_score = fields[2]
    	   taxa = hits.group(1)
    	   out_line = read + ", " + taxa + ", " + last_score + "\n"
    	   output_handle.write(out_line)

    output_handle.close()

    exit()

# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
