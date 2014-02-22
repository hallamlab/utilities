#!python

"""
fasta_remapper.py

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
            # print fasta_filename
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


    def create_database(in_fasta, exe, out_folder):
        name = in_fasta.split("/")[-1]
        name = re.sub(r'\..*', '',name)
        out_name = out_folder + "/" + name
        cmd='%s -dbtype nucl -in %s -out %s' %(exe, in_fasta, out_name)
        result = getstatusoutput(cmd)
        return out_name

    def blast(exe, query, database, percent_id, output_folder):
        query_name =  re.sub(r'\..*', '',query.split("/")[-1]) 
        db_name =  re.sub(r'\..*', '',database.split("/")[-1]) 
        out_name = output_folder + "/" + query_name + "_" + db_name + ".blastout"
        cmd = '%s -query %s -db %s -outfmt 6 -perc_identity %s -min_raw_gapped_score %s -out %s' %(exe, query, database, str(percent_id), str(100), out_name)
        # print cmd
        result= getstatusoutput(cmd)
        return out_name

    # helper function to get directories from a given path
    def get_directories(path):
        ret = []
        for f in os.listdir(path):
            if not os.path.isfile(os.path.join(path,f)):
                ret.append(f.strip())
        return ret

    def get_files(path):
        ret = []
        for f in os.listdir(path):
            if os.path.isfile(os.path.join(path,f)):
                ret.append(f.strip())
        return ret

    def parse_naming_file(file):
        try:
            handle = open(file, 'r')
        except IOError:
            print "Could not open mapping file."

        lines = handle.readlines()
        library_name_to_full_name = {}
        library_name_to_external = {}
        mappings = {"library_name_to_full_name":library_name_to_full_name, "library_name_to_external":library_name_to_external }
        lines = lines[2:len(lines)] # drop the header lines
        for l in lines:
            if re.match('.*\|.*\|.*\|.*\|.*', l):
                fields = l.split("|")
                lib_name = fields[1].strip()
                lib_name_full =  fields[2].strip()
                external_id = fields[3].strip()
                mappings["library_name_to_full_name"][lib_name] = lib_name_full
                mappings["library_name_to_external"][lib_name] = external_id

        return mappings



what_i_do = "Strip fasta file names"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_fasta', type=str, nargs='?',
                required=True, help='a selection of one or many fosmid input .fasta files (required)', default=None)
parser.add_argument('-m', dest='name_map', type=str, nargs='?',
            required=True, help='map file for the fasta headers (required)', default=None)                             
parser.add_argument('-o', dest='output_directory', type=str, nargs='?',
                required=True, help='the name a target output directory (required)', default=None)

def main(argv):
    args = vars(parser.parse_args())
    
    f = args["input_fasta"]
    reader = FastaReader(f)
    map_handle = open(args["name_map"], "r")
    map_lines = map_handle.readlines()
    map_dict = {}
    for l in map_lines:
        fields = l.split("\t") 
        map_dict[">" + fields[0].strip()] = ">" + fields[1]
    out_handle = open(os.path.abspath(args["output_directory"]) + "/" + os.path.basename(f) + ".new.fasta", "w")
    for record in reader:
        new_name = ""
        if record.name in map_dict:
            new_name = map_dict[record.name]
        else:
            new_name = record.name
        out_handle.write(new_name)
        out_handle.write(record.sequence + "\n")
    out_handle.close()
    exit()


# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
