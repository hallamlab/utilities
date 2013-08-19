#!/usr/bin/python
"""
do_something_script.py

Created by Niels Hanson on 2013-08-16.
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


what_i_do = "I'm a script that does something"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_files', type=str, nargs='+',
                required=True, help='a selection of one or many input .fasta files (required)', default=None)
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
               required=True, help='the name of the output frequeny matrix (required)', default=None)

# helper function to get directories from a given path
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
            close(file)
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
    
# Functional classes    
"""
A class to represent an individual .fasta record within a fasta file.
It contains three fields: (1) longname - representing the whole annotation 
contained after the carat symbol (>), (2) sequence - the amino acid or nucleotide 
sequence, and (3) name a shortened name based on the information at appears before
first space in the annotation line.
"""
class FastaRecord():
    def __init__(self, longname, sequence):
      self.longname = longname
      self.sequence = sequence
      fields = [ x.strip() for x in self.longname.split(' ') ]
      if len(fields) > 0:
         self.name = fields[0]
      else:
         self.name = None

"""
Takes a fasta file as input and returns a list of FastaRecords for
each correctly parsed record in the file.
"""
class FastaReader():
    stop = False
    START_PATTERN = re.compile(r'^>')
    name = None
    future_name =None
    sequence=""
    def __init__(self, fasta_filename):
        try:
            self.file = open(fasta_filename, 'r')
            print fasta_filename
        except IOError:
            print "Cannot open fasta file " + fasta_filename

    def __iter__(self):
        return self


    def next(self):
        if self.stop:
          raise StopIteration

        try:
           if not self.name: 
               self.name = self.file.readline().strip()
           line = self.file.readline().strip()
        except:
           line = None


        if not line:
           self.stop = True
           raise StopIteration


        fragments = []
        while line and not self.START_PATTERN.search(line):
            fragments.append(line.strip()) 
            line = self.file.readline()

       # print line
        if self.future_name:
            self.name = self.future_name

        if line:
          self.future_name = line.strip()

        self.sequence =''.join(fragments)
        self.seqname = self.name

        return FastaRecord(self.name, self.sequence)

"""
A function to calculate the frequency of valid tetramers in a given sequence. Returns a freq_hash
containing the tetramers and their frequencies. If given a frequency hash, frequenies found in the
sequence are added to the existing hash.
"""
def calculate_tetramer_freq(sequence, freq_hash, non_nucleotides, n_pattern):
    if sequence == None:
        print "what are you doing?"
        exit()
    if freq_hash == None:
        freq_hash = {}
    
    sequence = sequence.upper()
    sequence = non_nucleotides.sub("",sequence)
    length = len(sequence)
    
    for i in range(0,length-3):
        mer = sequence[i:i+4]
        if mer not in freq_hash:
            freq_hash[mer] = 0
        freq_hash[mer] = freq_hash[mer] + 1
        
    return freq_hash

# the main function
def main(): 
    # parse arguments
    args = vars(parser.parse_args())
    
    # setup our regular expresions we are going to use a lot
    non_nucleotides = re.compile("[^ATCG]")
    n_pattern = re.compile("N")
    
    # hash to store 4mer results from each input file
    files_to_4mers = {}
    
    # create a list of all possible 4mers
    nuc = ["A", "T", "C", "G"]
    mer_list = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    mer = nuc[i] + nuc[j] + nuc[k] + nuc[l]
                    mer_list.append(mer)
    
    # for each input fasta file extract sequences using the FastaReader and use calculate_tetramer_freq 
    # to sum the frequency of 4mers in each sequence. Storing the result in files_to_4mers hash
    for f in args["input_files"]:
        reader = FastaReader(f)
        freq_hash = None
        for record in reader:
            freq_hash = calculate_tetramer_freq(record.sequence, freq_hash, non_nucleotides, n_pattern)
        
        filename = f.split("/")[-1] # python tricks to get filename
        files_to_4mers[filename] = freq_hash # store calculated frequency hash
    
    # write out the results
    output_handle = open(args["output_file"], "w")
    # construct the header
    header = ""
    for sample in files_to_4mers:
        header = header + "\t" + sample
    output_handle.write(header + "\n")
    # construct the mer count accross samples
    line = ""
    for mer in mer_list:
        line = mer
        for f in files_to_4mers:
            num = 0
            if mer in files_to_4mers[f]:
                num = files_to_4mers[f][mer]
            line = line + "\t" + str(num)
        output_handle.write(line + "\n")
    
    # close file
    output_handle.close()
    exit()
    
    
# call the main function
if __name__ == "__main__":
    sys.exit(main())

