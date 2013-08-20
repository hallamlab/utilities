#!/usr/bin/python
"""

Created by Kan Cheung on 2013-07-04.
Copyright (c) 2013 Steven J. Hallam Laboratory. All rights reserved.
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


what_i_do = "This script takes in a list of edges (.csv) and create an upper triangular matrix."
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_files', type=str, nargs='+',
                required=True, help='a selection of input files (required)', default=None)
parser.add_argument('-o', dest='output_files', type=str, nargs='?',
               required=True, help='the output file (required)', default=None)

# helper function to get directories from a given path
def get_directories(path):
    ret = []
    for f in os.listdir(path):
        if not os.path.isfile(os.path.join(path,f)):
            ret.append(f.strip())
    return ret

def process_input_file(input_file):
    # lets do our checks
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
    
        myinput.close() # always close your files
        print "closed file " + input_file
        

def number_of_taxa(input):
    file_handle = open(input, 'r')
    header = file_handle.readline() 
    lines = file_handle.readlines()
    taxa = [] 
    for line in lines:
        line_array = line.split(",") 
        for element in line_array: 
            t = element.replace("\n", "") 
            if t not in taxa: 
                taxa.append(t.strip())
    return taxa

"""
Check the arguments given
"""
def check_arguments(arguments):
    for file in arguments["input_files"]:
        try:
            open(file)
            close(file)
        except IOError:
            print "Could not open " + file
            sys.exit(1)

def count_data_edges(input_file, dict_of_dict):
    if (input_file == None) or ( not isinstance(input_file, str)) or len(input_file) == 0:
        print "some problems with the arguments"
        return -1
    else:
        # okay there's something there we'll try to open it
        try:
            #print "open " + input_file
            file_handle = open(input_file,'r')
        except IOError:
            print "Cannot open " + str(input_file)
            exit()
            
        header = file_handle.readline()
        lines = file_handle.readlines()
        for line in lines:
            line_array = line.split(",")
            # add to the dictionary of dictionaries
            taxa1 = line_array[0].replace("\n","").strip()
            taxa2 = line_array[1].strip().replace("\n","")
            dict_of_dict[taxa1][taxa2] = dict_of_dict[taxa1][taxa2] + 1
        return dict_of_dict
    
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
    
# the main function
def main(): 
    # parse arguments
    args = vars(parser.parse_args())

    # calculate N
    taxa_list = number_of_taxa(args["input_files"][0])
    
    # dictionary of dictionaries
    dict_of_dict = {}
    for t1 in taxa_list:
        dict_of_dict[t1] = {} # prints out the relationships t1 has with all other taxa
        for t2 in taxa_list:
            dict_of_dict[t1][t2] = 0

    # fill dictionary with taxa
    counted_dictionary = count_data_edges(args["input_files"][0], dict_of_dict)
    
    i = 0
    j = 0
    N = len(taxa_list)

    # add lower trianguar to upper triangular
    for i in range(N): 
        for j in range(i+1,N): 
            counted_dictionary[taxa_list[i]][taxa_list[j]] = counted_dictionary[taxa_list[i]][taxa_list[j]] + counted_dictionary[taxa_list[j]][taxa_list[i]] 
            
    # write out our results as a csv file        
    write_to_csv(counted_dictionary, taxa_list, "output_file.csv")

    
# call the main function
if __name__ == "__main__":
    sys.exit(main())


    
    
    
                

