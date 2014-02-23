#!/usr/bin/python

"""
flipped.py

Created by Mike Wu
Edited by Niels Hanson
"""

try:
   import csv
   import argparse
   import os
   import re
   import sys
except:
   print "Failed to load some libraries"
   exit(1)

what_i_do = "I create a non-redundant list of edges from a Network Edge file."
parser = argparse.ArgumentParser(description=what_i_do)

# add arguments to the parser
parser.add_argument('-e', dest='edge_file', type=str, nargs='?',
                required=True, help='input edge file', default=None)
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
                required=True, help='the reulting output file', default=None)


def main(argv):
   args = vars(parser.parse_args())
   seen_list = []
   
   file = open(args["edge_file"], "r+")
   out = open(args["output_file"], "w")

   f = csv.reader(file, delimiter='\t') 
   o = csv.writer(out, delimiter='\t') 

   def flipString(s):
	return [s[1] , s[0]]

   for row in f:
	if flipString(row) in seen_list:
	   row = flipString(row)
	else:
	   seen_list.append(row)
	o.writerow(row)

if __name__ == "__main__":
   main(sys.argv[1:])
