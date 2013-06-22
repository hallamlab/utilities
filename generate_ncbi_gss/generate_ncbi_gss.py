#!/usr/bin/python
# File created May 2013
from __future__ import division

__author__ = "Niels W. Hanson"
__copyright__ = "Copyright 2013"
__credits__ = ["r"]
__version__ = "1.0"
__status__ = "Release"

try:
     import re
     import sys
     import os
     import pdb
     from optparse import OptionParser, OptionGroup
except:
     print """ Could not load some user defined module functions"""
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

# Script Usage Info
script_name = "generate_ncbi_gss.py"
usage= script_name + """ --fasta_dir fasta_directory [-p] gss_parameter_file -o output.gss """
parser = OptionParser(usage)

# Possible Script Arguments
parser.add_option( "--fasta_dir", dest="fasta_dir", help='folder containing a number of .fasta files')
parser.add_option( "-o", dest="output_file", help='the output .gss file')
parser.add_option( "-p", dest="parameter_file", help="a parameter file for .plc file creation")

"""Check the arguments and proceeds if requirements met"""
def check_arguments(opts, args):
    if opts.fasta_dir == None:
         print """REQUIRED: You must provide a directory full of .fasta files"""
         return False
         
    if opts.output_file == None:
         print """REQUIRED: You must have an output file."""
         return False
         
    if opts.parameter_file == None:
         print """Note: Provide a parameter file for a .plc companion file."""
         
    return True

"""returns a list of (non-hidden) files from a given path"""
def get_files(path):
    ret = []
    for f in os.listdir(path):
        if os.path.isfile(os.path.join(path,f)):
            if not re.match("^\.", f.strip()):
                ret.append(f.strip())
    return ret

"""returns a list of directories from a given path"""
def get_directories(path):
    ret = []
    for f in os.listdir(path):
        if not os.path.isfile(os.path.join(path,f)):
            ret.append(f.strip())
    return ret

"""
Performs a variety of quality control steps on sequences.
"""
def qc_sequence(input_seq, N=100, L_min=50, L_max=900, l=2, m=10, n=2, p=5):
    
    # # trim ambiguous bases from head and tail N characters
    input_seq = trim_ambig_head_tail(input_seq, N)
    # # check if trimmed sequence in acceptable range
    if (len(input_seq) > L_max or len(input_seq) < L_min):
        print "Trimmed sequence not in accepted range: (L_min= " + str(L_min) + ", L_max= " + str(L_max) + ")"
        return None
    # # check remaining characters for ambiguous characters
    if violate_percent_ambig(input_seq, p):
        print "Trimmed sequence has too many ambiguous characters (not ATCG): (p= " + str(p) + ")"
        return None
    # # detects for acceptable amount of heterpolymers
    if violate_hetero_polymers(input_seq, m, n, l):
        print "Sequence has too many hetero_polymers: (l, m, n) = (" + str(l) + ", " + str(m) + ", " + str(n) + " )"
        return None
    # sequence passed all qc checks
    return input_seq

"""
Trims the first and last N characters for ambiguous nucleotides (not ATCG)
"""
def trim_ambig_head_tail(seq, N):
    seq_len = len(seq)
    if seq_len <= N:
        print "Sequence too short to trim"
        return seq
    left_N = seq[0:N]
    right_N = seq[(seq_len-N):seq_len]
    
    non_nucleotides = re.compile("[^ATCGatcg]+")
    l_results = non_nucleotides.finditer(left_N)
    l_max_index = 0
    if l_results:
        for l_match in l_results:
            if l_match.end(0) > l_max_index:
                l_max_index = l_match.end(0)
    r_results = non_nucleotides.finditer(right_N)
    r_min_index = seq_len
    if r_results:
        for r_match in r_results:
            if r_match.start(0) < r_min_index:
                r_min_index = r_match.start(0)
    return seq[l_max_index:seq_len-N+r_min_index]

"""
Checks to see if percentage of ambiguous characters with within percent theshold p
"""
def violate_percent_ambig(seq, p):
    non_nucleotides = re.compile("[^ATCGatcg]")
    matches = non_nucleotides.findall(seq)
    percent_ambig = len(matches)/len(seq) * 100
    if percent_ambig > p:
        return True
    else:
        return False

"""
Checks to see if there are too many heteropolymer sets based on three parameters:
l = the length of heteropolymer
m = the maximum number of repeats of a polymer per violation
n = the maximum number of violations per sequence
"""
n_curr = 0 # global variable

def violate_hetero_polymers(seq, m, n, l):
    repeat_hash = {} # keeps track of length l repeats and rangers
    global n_curr
    n_curr = 0 # number of current violations
    critical = False # flag to reject sequence if n_curr > n
    
    # iterate through sequence
    for i in range(len(seq)):
        s = i - l
        if s < 0:
            s = 0
        critical = update_keys(seq[s:i], repeat_hash, i, n, m)
        if critical == True:
            return True
    
    # one final pass to catch non-broken repeats near end
    for repeat in repeat_hash:
        arr = repeat_hash[repeat].split(":")
        if check_violation(int(arr[0]), int(arr[1]), repeat, m):
            n_curr += 1
            if n_curr > n:
                # critical violation
                return True
        
    
    return False


"""
Given a subsequence of a nucleotide string, sub_seq, function determines based on
contents of the repeat_hash current number of violocations, n_curr, and 
size of violation m, if a critical violation occurs in the current substring.
"""
def update_keys(sub_seq, repeat_hash, curr_pos, n, m):
    global n_curr
    sub_len = len(sub_seq)
    if sub_len == 0:
        return False
    for i in xrange(sub_len, -1 ,-1):
        sub_window = sub_seq[i:sub_len]
        if len(sub_window) == 0:
            continue
        if len(sub_window) > 1 and len(set(sub_window)) == 1:
            # string window only contains one character
            # extend individual character in the hash
            sub_window = iter(set(sub_window)).next()
            curr_range = repeat_hash[sub_window]
        else:
            if sub_window not in repeat_hash:
                repeat_hash[sub_window] = str(curr_pos - (sub_len - i)) + ":" + str(curr_pos)
            curr_range = repeat_hash[sub_window]
        arr = curr_range.split(":")
        j = int(arr[0])
        k = int(arr[1])
        l = curr_pos - (sub_len - i)
        p = curr_pos
        if (set(range(j,k) + range(l,p)) == set(range(j,p))):
            # range is contiguous
            repeat_hash[sub_window] = str(j) + ":" + str(p)
        else:
            if check_violation(j,k,sub_window, m):
                n_curr += 1
                if n_curr > n:
                    # critical violation
                    return True
            repeat_hash[sub_window] = str(l) + ":" + str(p)
        
    return False    
"""
Does a quick calculation to see if there is a violation in the current
sub-window. Returning True if so and False otherwise
"""
def check_violation(i,j,key,m):
    if ((j-i) / len(key)) > m:
        return True
    else:
        return False


"""
Handles a list of fasta files and set of plc instructions (optional) 
and produces a GSS file for each.
"""
def handle_fasta_files(fasta_files,opts,plc_options=None):
    for fasta_file in fasta_files:
        # clean and prepare variables
        cur_dur = opts.fasta_dir + "/" + fasta_file
        library = re.sub("\.(fa|fas|fasta|fna|f)$","", fasta_file)
        
        # open the output file (one per input fasta)
        try:
            myout = open(library + ".gss",'w')
        except IOError:
            print "Cannot open " + str(output_file)
        
        # open the current file
        reader = FastaReader(cur_dur)
        
        # for each FastaRecord in file
        for record in reader:
            # pull the short_name and sequence
            seq_id = record.name
            seq_id = seq_id.replace(">","")
            seq = record.sequence
            seq = qc_sequence(seq)
            if seq != None:
                # sequence passed quality control
                # generate the GSS record 
                if plc_options != None:
                    line = ""
                    line += "TYPE: GSS\n"
                    line += "STATUS: New\n"
                    line += "CONT_NAME: " + plc_options["CONT_NAME"] + "\n"
                    line += "CITATION: " + "\n" + plc_options["PUB_TITLE"] +"\n"
                    line += "LIBRARY: " + library + "\n"
                    line += "GSS#: " + seq_id + "\n"
                    line += "CLASS: fosmid ends" + "\n"
                    line += "SEQUENCE: " + seq + "\n"
                    line += "\n"
                    line += "||\n\n"
                else:
                    line = ""
                    line += "TYPE: GSS\n"
                    line += "STATUS: New\n"
                    line += "CONT_NAME: Hallam, S.J.\n"
                    line += "CITATION: " + "\n" + "A potential role for Marine Group A bacteria in the marine sulfur cycle" +"\n"
                    line += "LIBRARY: " + library + "\n"
                    line += "GSS#: " + seq_id + "\n"
                    line += "CLASS: fosmid ends" + "\n"
                    line += "SEQUENCE: " + seq + "\n"
                    line += "\n"
                    line += "||\n\n"
                
                # write out the current record
                myout.write(line)
            
        # close output gss file
        myout.close()

"""
Create PLC option file
"""    
def generate_plc_file(fasta_files, plc_option_file):
    plc_options = {}
    # open the output file (one per input fasta)
    try:
        myplc_file = open(plc_option_file,'r')
    except IOError:
        print "Cannot open " + str(output_file)
        exit()
    
    lines = myplc_file.readlines()
    myplc_file.close()
    plc_pattern = re.compile(r'((PUB.*?|LIB.*?|CONT.+?).*?):(.*$)')
    
    for line in lines:
        hits = plc_pattern.search(line)
        if hits:
            plc_options[hits.group(1).strip()] = hits.group(3).strip()
            
    plc_output = open(plc_option_file + ".plc.txt", "w")
    # create plc file
    # Publication (.pub)
    pub = "TYPE: Pub" + "\n"
    pub += "TITLE: " + plc_options["PUB_TITLE"] + "\n"
    pub += "AUTHORS: " + plc_options["PUB_AUTHORS"] + "\n"
    pub += "JOURNAL: " + plc_options["PUB_JOURNAL"] + "\n"
    pub += "VOLUME: " + plc_options["PUB_VOLUME"] + "\n"
    pub += "ISSUE: " + plc_options["PUB_ISSUE"] + "\n"
    pub += "PAGES: " + plc_options["PUB_PAGES"] + "\n"
    pub += "YEAR: " + plc_options["PUB_YEAR"] + "\n"
    pub += "STATUS: " + plc_options["PUB_STATUS"] + "\n"
    pub += "||" + "\n" + "\n"
    
    plc_output.write(pub)
    
    # Library (.lib)
    for fasta_file in fasta_files:
        lib = ""
        library = re.sub("\.(fa|fas|fasta|fna|f)$","", fasta_file)
        lib += "TYPE: Lib" + "\n"
        lib += "NAME: " + library + "\n"
        lib += "VECTOR: " + plc_options["LIB_VECTOR"] + "\n"
        lib += "V_TYPE: " + plc_options["LIB_V_TYPE"] + "\n"
        lib += "DESCR: " + plc_options["LIB_DESCR"] + "\n" 
        lib += "||" + "\n" + "\n"
        
        plc_output.write(lib)
        
    # Contact (.cont)
    cont = "TYPE: Cont" + "\n"
    cont += "NAME: " + plc_options["CONT_NAME"] + "\n"
    cont += "TEL: " + plc_options["CONT_TEL"] + "\n"
    cont += "EMAIL: " + plc_options["CONT_EMAIL"] + "\n"
    cont += "LAB: " + plc_options["CONT_LAB"] + "\n"
    cont += "INST: " + plc_options["CONT_INST"] + "\n"
    cont += "ADDR: " + plc_options["CONT_ADDR"] + "\n"
    cont += "||" + "\n" + "\n"
    
    plc_output.write(cont)
    plc_output.close()
    return plc_options


"""the main function"""
def main(argv): 
    # grab and check command line arguments
    (opts, args) = parser.parse_args() 
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    
    fasta_files = [] # list of fasta_files
    
    # get rid of all non-fasta files by suffex
    for file in get_files(opts.fasta_dir):
        if re.search("\.(fa|fas|fasta|fna|f|faa)$", file):
            # file is fasta
            fasta_files.append(file)
        else:
            print "Warning " + file + " skipped file -- non-regular (.fa,.fas,.fasta,.fna,faa) suffix"
    
    # generate companion plc filed
    if opts.parameter_file != None:
        plc_options = generate_plc_file(fasta_files, opts.parameter_file)
    
    # handle fastafiles file and create gcc files
    handle_fasta_files(fasta_files,opts, plc_options)


# call the main function with command line argument vector
if __name__ == "__main__":
    main(sys.argv[1:])



