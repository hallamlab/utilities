#!python

"""
compare.py
Written by Mike Wu
Edited by Niels Hanson
"""


try:
   import csv
   import numpy as np
   import sys
   import argparse
except:
   print "Could not load some modules"

what_i_do = "Count the number of edges between taxa." 
parser = argparse.ArgumentParser(description=what_i_do)

# add arguments to the parser
parser.add_argument('-i', dest='input_flipped_file', type=str, nargs='?',
                required=True, help='input flipped edge file from flip.py (required)', default=None)
parser.add_argument('-s', dest='schema_file', type=str, nargs='?',
                required=True, help='a file specifying the order of the taxa in a column', default=None)             
parser.add_argument('-o', dest='output_matrix', type=str, nargs='?',
                required=True, help='output matrix of counted edges (required)', default=None)

counts = {}
def parseFile(f):
    file = open(f,"r")
    lines = file.readlines()

    for l in lines:
        l = l.strip()
        split_up = split(l) # create variable to hold split up string
        first = split_up[0] # first thing in the split
        second = split_up[1]

        if not counts.has_key(first):
            counts[first] = {}
        if not counts[first].has_key(second):
            counts[first][second] = 0

        counts[first][second] = counts[first][second]+1 

def split(str):
	return str.split("\t")

l = list()
def getKeys(args):
    schema = open(args["schema_file"],"r")
    lines = schema.readlines()

    for line in lines:
        l.append(line.strip())

rows = dict()
def generateRows():
    for k,v in counts.iteritems():
        if not rows.has_key(k):
        	rows[k] = dict()
        rows[k] = v

writtenData = list()
def output():
    for header in l:
        rowEntry = list()
        if header not in rows:
            rowEntry = [ 0.0 for x in range(len(l)) ]
            writtenData.append(rowEntry)
        else:
            for rk in rows:
                if rk == header:
                    for colHeaders in l:
                        if colHeaders in rows[rk]:
                            rowEntry.append(rows[rk][colHeaders])
                        else:
                            rowEntry.append(0)                
                    writtenData.append(rowEntry)

def upperTriangular(args):
	numRows = len(writtenData)-1
	numCols = numRows # assumes output is square

	for i in range(numRows):
		for j in range(numCols):
			if numRows-i>-1 and numCols-j>-1 and i!=j:
				# print 'at ', i, j, numRows, numCols,' adding ', writtenData[i][j], ' with ', writtenData[numRows-i][numCols-j]
 				writtenData[i][j] += writtenData[j][i]
	out = np.matrix(writtenData)
	out = np.triu(out)
	np.savetxt(args["output_matrix"], out, fmt='%i',delimiter="\t")


def main (argv):
   
   args = vars(parser.parse_args())
   parseFile(args["input_flipped_file"])
   getKeys(args)
   generateRows()
   output()
   upperTriangular(args)

if __name__ == "__main__":
   main(sys.argv[1:])

