#!python

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

what_i_do = """Prepare a GreenGenes fasta sequence database suitable for BLAST/LAST processing. 
Sequences and files can be obtained from http://greengenes.secondgenome.com/"""
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser

parser.add_argument('-f', dest='fasta_sequences', type=str, nargs='?',
          required=True, help='raw fasta sequences with GreenGene IDs in the headers', default=None)
parser.add_argument('-t', dest='taxa_map', type=str, nargs='?',
          required=True, help='taxonomy map file connecting GreeGene IDs to the taxonomy', default=None)
parser.add_argument('-a', dest='accession_map', type=str, nargs='?',
        required=False, help='taxonomy map file connecting GreeGene IDs to the taxonomy', default=None)
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
         required=True, help='name of the prepared database .fasta file', default=None)


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
Process ggid to taxa file
"""
def process_taxa(taxa_file):
    try:
        taxa_handle = open(taxa_file, "r")
    except:
        print "problem opening taxa file " + taxa_file
    lines = taxa_handle.readlines()
    taxa_handle.close()
    
    gg_to_taxa = {}
    for l in lines:
        fields = l.split("\t")
        gg_to_taxa[fields[0]] = fields[1].strip("\n")
    
    return gg_to_taxa


"""
process ggid to accession
"""
def process_accession(acc_file):
    try:
        accession_handle = open(acc_file, "r")
    except:
        print "problem opening accession map file " + acc_file
    lines = accession_handle.readlines()
    accession_handle.close()
    
    gg_to_accession = {}
    for l in lines:
        fields = l.split("\t")
        gg_to_accession[fields[0]] = fields[2].strip("\n")
    
    return gg_to_accession


# the main function
def main():
    # parse arguments
    args = vars(parser.parse_args())
    # read taxonomy
    
    gg_to_taxa = process_taxa(args["taxa_map"])
    if args["accession_map"]:
        gg_to_accession = process_accession(args["accession_map"])
    
    # read the input fasta file
    reader = FastaReader(args["fasta_sequences"])
    out_file = open(args["output_file"], "w")
    for record in reader:
        temp = record.name.lstrip(">")
        line = ">" + temp
        if temp not in gg_to_accession:
            print "not found accession"
        if temp not in gg_to_taxa:
            print "not found in taxa"
        line = line + " " + gg_to_accession[temp] + " " + gg_to_taxa[temp] + "\n"
        out_file.write(line)
        out_file.write(record.sequence + "\n")
        
    out_file.close()
    
    exit()
    
# call the main function
if __name__ == "__main__":
    sys.exit(main())

# >16 AF056938.1 Methanocaldococcus fervens str. AG86; DSM4213 k__Archaea; p__Euryarchaeota; c__Methanococci; o__Methanococcales; f__Methanocaldococcaceae; g__Methanocaldococcus;