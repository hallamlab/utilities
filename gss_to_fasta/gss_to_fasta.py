#!/usr/bin/python
"""
do_something_script.py

Created by Niels Hanson on 2013-04-21.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

# load some packages
from __future__ import division
try:
     import re
     import sys
     import os
     import types
     import argparse
except:
     import_exception = """ Could not load some modules """
     print import_exception 
     sys.exit(3)


what_i_do = """ This script converts NCBI Genome Survey Sequences files (.gss) files to .fasta files by 
extracting all the sequences
"""
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_files', type=str, nargs='+',
                required=True, help='a selection of input files (required)', default=None)
parser.add_argument('-o', dest='output_files', type=str, nargs='?',
               required=True, help='the output file (required)', default=None)
 
"""
Helper function to get directories from a given path
"""
def get_directories(path):
    ret = []
    for f in os.listdir(path):
        if not os.path.isfile(os.path.join(path,f)):
            ret.append(f.strip())
    return ret
    
def process_input_file(input_file):
    # lets do our checks
    # print "In process_input_file()"
    if (input_file == None) or ( not isinstance(input_file, str)) or len(input_file) == 0:
        print "some problems with the arguments"
        return -1
    else:
        # okay there's something there we'll try to open it
        try:
            print "open " + input_file
            myinput = open(input_file,'r')
        except IOError:
            print "Cannot open " + str(input_file)
            exit()
        
        header = myinput.readline().strip("\n").split(",")
        data_matrix = {}
        
        myinput.close() # always close your files
        print "closed file " + input_file

"""
Check the arguments given
"""
def check_arguments(arguments):
    for file in arguments["input_files"]:
        try:
            open(file)
        except IOError:
            print "Could not open " + file
            sys.exit(1)
    
def write_to_csv(counted_dictionary, taxa_list, output_file_name):

    # open an output file
    try:
        output_handle = open(output_file_name, "w")
    except IOError:
        print "Could not open " + output_file_name
        sys.exit(1)
        
    # create the header
    header = ""
    
    for taxa in taxa_list:
        header = header + ", " + taxa 
    header = header + "\n"
    output_handle.write(header)
    # print header
    
    i = 0
    j = 0
    N = len(taxa_list)
    
    line = ""
    for i in range(N):
        # add name to row
        line = taxa_list[i]
        for j in range(N):
            # upper tri
            if j >= i:
                line = line + ", " + str(counted_dictionary[taxa_list[i]][taxa_list[j]]) 
            else:
                line = line + ", 0"
                # lower tri add a zero
        line = line + "\n"
        output_handle.write(line)
    
    # close output file
    output_handle.close()
    
def gss_to_fasta(input_files):
    GSS_header_pattern = re.compile("^GSS#:(.*)\n") # ^: beginning of line || . : any char || * : any number of previous pattern || (): pattern capture
    sequence_start = re.compile("^SEQUENCE:")
    sequence_stop = re.compile("^\|\|")
    
    # run through the files
    for i in input_files:
        # open the file, read lines, and close
        handle = open(i, "r")
        lines = handle.readlines()
        handle.close()
        
        # open output file
        out_handle = open(i + ".fasta", "w")
        
        # header for each sequence and variables for current sequence
        fasta_header = ""
        current_sequence = []
        collect = False
        for j in lines:
            hits = GSS_header_pattern.search(j)
            if hits:
                fasta_header = hits.group(1).strip()
                line = ">" + fasta_header + "\n"
                out_handle.write(line)
            if sequence_start.match(j):
                collect = True
            if sequence_stop.match(j):
                collect = False
                current_sequence = "".join(current_sequence)
                current_sequence = re.sub("SEQUENCE: ", "",current_sequence)
                current_sequence = current_sequence.strip()
                current_sequence = current_sequence + "\n"
                out_handle.write(current_sequence)
                current_sequence = []
            if collect:
                current_sequence.append(j)
                
        out_handle.close()        
                
            
        
                
                       


# the main function
def main(): 
    # parse arguments
    args = vars(parser.parse_args())
    gss_to_fasta(args["input_files"])
    
    

    
    
# call the main function
if __name__ == "__main__":
    sys.exit(main())


    
    
    
                

