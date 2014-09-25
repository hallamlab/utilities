#!/usr/bin/python
# File created on Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar, Niels W. Hanson"
__copyright__ = "Copyright 2013, MetaPathways"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import re
     import sys
     import os
     from optparse import OptionParser, OptionGroup
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwayrc\""""
     print """ """
     sys.exit(3)

script_name = "cazy_foly_summary.py"
usage= script_name + """ --parse_blast blast_directory [-n 1] -o output.csv"""
parser = OptionParser(usage)

# two arguments the parsed-blast outputs and the location to place the output file
parser.add_option( "--parse_blast", dest="parse_blast", help='folder containing the parse BLAST files for CAZy and FOLy')
parser.add_option("-n", "--num_hits", dest="num_hits", default=1, help="number of hits per query sequence to count")
parser.add_option( "-o", dest="output_file", help='the output file')


# check the arguments
def check_arguments(opts, args):
    
    if opts.parse_blast == None:
         print """Must have the \"parse_blast\" directory"""
         return False
         
    if opts.output_file == None:
         print """Must have an output file"""
         return False
         
    return True

# get files from a given path
def get_files(path):
    ret = []
    for f in os.listdir(path):
        if os.path.isfile(os.path.join(path,f)):
            if not re.match("^\.", f.strip()):
                ret.append(f.strip())
    return ret

# get directories from a given path
def get_directories(path):
    ret = []
    for f in os.listdir(path):
        if not os.path.isfile(os.path.join(path,f)):
            ret.append(f.strip())
    return ret

# parse a cazy blast output for families
def parse_cazy(cazy_file, num_hits):
    # open the file
    try:
        my_file = open(cazy_file,'r')
    except IOError:
        print "Cannot open " + str(my_file)
        
    # read in the lines and close
    lines=my_file.readlines()
    my_file.close()
    
    # cazy family dictionary
    cazy_families = {}
    query_count = {} # keeps track of read counts
    cazy_pattern = re.compile("\(([A-z]+)\:([0-9_]+)\)")
    for line in lines:
        words = [x.strip() for x in line.rstrip().split("\t") ]
        hits = cazy_pattern.search(words[8])
        if hits:
            cur_read = str(words[0].strip())
            cur_fam = str(hits.group(1)) + str(hits.group(2))
            if cur_read not in query_count:
                query_count[cur_read] = 0
            query_count[cur_read] += 1
            if query_count[cur_read] <= num_hits:
                if cur_fam not in cazy_families:
                    cazy_families[cur_fam] = 0
                cazy_families[cur_fam] += 1
    
    return cazy_families

# given the location of the BLAST directory it will find all parsed output files with
# CAZy in the file name
def parse_cazy_families(files, num_hits):
    
    # isolate parsed CAZy files from the BLAST directory
    cazy_file_list = {}
    for f in files:
        if ( re.match(".*cazy.*", f, re.IGNORECASE) and re.match(".*blastout\.parsed.*", f, re.IGNORECASE)):
            cazy_file_list[d].append(d + blast_dir + f)
    
    # parse families for each sample
    cazy_results = {}

    for s in cazy_file_list:
        for f in cazy_file_list[s]:
            target_dir = base_dir + "/" + f
            cazy_results[s] = parse_cazy(target_dir, num_hits)
            
    return cazy_results

# parse a foly blast output for families
def parse_foly(foly_file, num_hits):
    # open the file
    try:
        my_file = open(foly_file,'r')
    except IOError:
        print "Cannot open " + str(my_file)

    # read in the lines and close
    lines = my_file.readlines()
    my_file.close()

    # FOLy family dictionary
    foly_families = {}
    query_count = {} # keeps track of read counts
    foly_pattern = re.compile("\(\((.*)_(.*)\)\)")
    for line in lines:
        words = [x.strip() for x in line.rstrip().split("\t") ]
        hits = foly_pattern.search(words[8])
        if hits:
            cur_read = str(words.strip())
            cur_fam = str(hits.group(1))
            if cur_read not in query_count:
                query_count[cur_read] = 0
            query_count[cur_read] += 1
            if query_count[cur_read] <= num_hits:
                if cur_fam not in foly_families:
                    foly_families[cur_fam] = 0
                foly_families[cur_fam] += 1

    return foly_families


def parse_foly_families(base_dir, num_hits):
    
    # isolate parsed FOLy files from BLAST directory
    blast_dir = "/blast_results/"
    dir_list = get_directories(base_dir)
    foly_file_list = {}
    for d in dir_list:
        foly_file_list[d] = []
        file_list = get_files(base_dir + "/" + d + "/" + blast_dir)
        # find all the parsed lists
        
        for f in file_list:
            if ( re.match(".*foly.*", f, re.IGNORECASE) and re.match(".*blastout\.parsed.*", f, re.IGNORECASE)):
                foly_file_list[d].append(d + blast_dir + f)
                
    # parse families for each sample
    foly_results = {}

    for s in foly_file_list:
        for f in foly_file_list[s]:
            target_dir = base_dir + "/" + f
            foly_results[s] = parse_foly(target_dir, num_hits)

    return foly_results    

def write_results(cazy_results, foly_results, output_file):
    # try to open the output file
    try:
        myout = open(output_file,'w')
    except IOError:
        print "Cannot open " + str(output_file)
    
    samples = []
    families = []
    # prepare ordinate
    for s in cazy_results:
        if s not in samples:
            samples.append(s)
        for f in cazy_results[s]:
            if f not in families:
                families.append(f)
    
    for s in foly_results:
        if s not in samples:
            samples.append(s)
        for f in foly_results[s]:
            if f not in families:
                families.append(f)
                
                
    samples.sort()
    families.sort()
    
    result_matrix = {}
    for i in samples:
        result_matrix[i] = {}
        for f in families:
            result_matrix[i][f] = {}
            if f in cazy_results[i]:
                result_matrix[i][f] = cazy_results[i][f]
            else:
                result_matrix[i][f] = 0
    
    # print formatted results
    header = ""
    for i in samples:
        header = header + "," + str(i)
    myout.write(header + "\n")
    print header
    
    for j in families:
        line = str(j)
        for i in samples:
            line = str(line) + "," + str(result_matrix[i][j])
        myout.write(line + "\n")
        print line
    
    myout.close()

# the main function
def main(argv): 
    # grab and check arguments
    (opts, args) = parser.parse_args()
    if not check_arguments(opts, args):
       print usage
       sys.exit(0)
    
    num_hits = int(opts.num_hits)
    
    # parse cazy and foly family counts from the blast output
    cazy_results = parse_cazy_families(opts.parse_blast, num_hits)
    foly_results = parse_foly_families(opts.parse_blast, num_hits)
    
    # print out results
    write_results(cazy_results,foly_results,opts.output_file)

# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

