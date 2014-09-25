#!/usr/bin/python
"""
fosmid_qc.py

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
     from sys import path
     from libs.metapaths_utils  import parse_command_line_parameters, fprintf
     from libs.sysutil import getstatusoutput, pathDelim
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)

what_i_do = "Fosmid Quality Control. Removes vector sequences and annotates with respect to available fosmid ends."
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_fastas', type=str, nargs='+',
                required=True, help='a selection of one or many fosmid input .fasta files (required)', default=None)             
parser.add_argument('-o', dest='output_directory', type=str, nargs='?',
                required=True, help='the name a target output directory (required)', default=None)
parser.add_argument('-d', dest='remove_me_fasta', type=str, nargs='?',
                required=True, help='a fasta file containing sequences to find and remove from input (required)', default=None)
parser.add_argument('-e', dest='fosmid_ends', type=str, nargs='+',
                required=True, help='location of the directory of fosmid-end fasta files for naming (required)', default=None)
parser.add_argument('-n', dest='naming_file', type=str, nargs='?',
                required=False, help='name mapping file linking Library_Name, Library FullName, and External_Identifier', default=None)
parser.add_argument('-p', dest='percent_id', type=str, nargs='?',
                required=False, help='minimum percent identity cutoff for blast homology searches (default 98)', default='98')
parser.add_argument('-l', dest='min_length', type=str, nargs='?',
                required=False, help='minimum length for sequences (default 2000)', default='2000')                
parser.add_argument('-b', dest='blast_executable', type=str, nargs='?',
                required=False, help='location of the blastn executable, will assume in PATH if not specified', default='blastn')
parser.add_argument('-f', dest='database_executable', type=str, nargs='?',
                required=False, help='location of the makeblastdb executable, will assume in PATH if not specified', default='makeblastdb')
parser.add_argument('--qc', dest='just_qc', action='store_true',
                required=False, help='flag to exit after the QC step', default=False)
               
def check_arguments(args):
   if args['input_fastas'] == None or args['output_directory'] == None or args['remove_me_fasta'] == None:
      # TODO check to see if these files/locations are actually accessable
      return True
   else:
      return False

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


# global convenience variables
input_fastas = None # list of input fasta filenames
input_dir = None # input fasta directory
remove_me_fasta = None # input remove_me fasta filename
remove_me_dir = None # directory
blast_exe = None # path to blastn
blast_db_exe = None # path to makeblastdb
naming_file = None # name of GSC mapping file
naming_file_dir = None # directory of GSC mapping file
fosmid_ends = None # list of fosmid end files
fosmind_end_dir = None # fosmid end directory
output_dir = None # output directory
temp_blastdb_dir = "temp_blastdb" # place to put blast databases
blastdb_list = [] # names of blastdbs
temp_blastout_dir = "temp_blastout" # place to put blast databases
blastouts = [] # names of blastouts
temp_fasta_dir = "temp_fasta"
temp_fastas = [] # names of temp_fastas
sep = "/" # compatability for windows users
fosmid_lengths = {} # lengths of all fosmids
xout_lengths = {} # lengths of all processed xout fosmids


def remove_primer_sequences(fasta_to_blastout, length_cutoff = 2000):
    for fasta in fasta_to_blastout:
        try:
            blast_handle = open(fasta_to_blastout[fasta], 'r')
        except IOError:
            print "Cannot open blastout file " + fasta_to_blastout[fasta]
            exit()
        
        lines = blast_handle.readlines()
        blast_handle.close()
        # "2014613006	pCC1FOS	100.00	953	0	0	464	1416	126	1078	0.0	1760"
        # pull out the hits
        targets = {}
        for line in lines:
            splits = line.split("\t")
            q_seq = splits[0].strip()
            if q_seq not in targets:
                targets[q_seq] = []
            q_start = splits[6].strip()
            q_end = splits[7].strip()
            coords = None
            if q_start > q_end:
                coords = [q_end,q_start]
            else:
                coords = [q_start, q_end]
            targets[q_seq].append(coords)
        # find targets in fasta
        reader = FastaReader(input_dir + sep + fasta)
        # open output_fasta
        try:
            fasta_handle = open(temp_fasta_dir + sep + fasta + ".xout", "w")
        except:
            "Problem opening fasta file"
        for record in reader:
            name = record.name.replace(">", "")
            seq = record.sequence
            if name in targets:
                qc_seq = seq
                for t in range(len(targets[name])):
                    i = int(targets[name][t][0])-1
                    j = int(targets[name][t][1])
                    l = j-i
                    if l > 20 and l < fosmid_lengths[fasta][name]:
                        # non-trival hit
                        # replace with X's
                        qc_seq = qc_seq[:i] + 'N'*(j-i) + qc_seq[j:]
                        
                        
                
                seq = qc_seq.strip('N') # remove end hits
            if len(seq) > length_cutoff:
                fasta_handle.write(">" + name + "\n")
                fasta_handle.write(seq + "\n")
        
        fasta_handle.close()


"""
Takes a list of filepaths and returns filenames
"""
def paths_to_filenames(in_list):
    for i in range(len(in_list)):
        in_list[i] = in_list[i].split(sep)[-1]
    return in_list


"""
From a list of fasta files, creates a hash of the lengths of
all sequences inside.
"""
def fasta_lengths(fasta):
    length_hash = {}
    reader = FastaReader(fasta)
    for r in reader:
        name = r.name.replace(">", "")
        l = len(r.sequence)
        if name not in length_hash:
            length_hash[name] = l
    return length_hash


# f_end = re.sub('\.fasta|\.fas|.faa|\.fa|\.fna|\.f','', f_end)
"""
Takes and evaluates a directory path, collects a list of files and removes them using
os.remove()
"""
def clean_dir(dir):
    for f in get_files(dir):
        path = dir + sep + f
        os.remove(path)

"""
process fosmid_hits for significant results
"""
def identify_fosmid_matches(foz_to_end_blastout):
    
    tail_pattern = re.compile(".*(_.+)")
    dot_pattern = re.compile("(^.+?)\.(.+)_")
    fb_pattern = re.compile(".*F.*", re.IGNORECASE)
    
    final_mapping = {}
    
    for f in foz_to_end_blastout:
        for j in foz_to_end_blastout[f]:
            try:
                handle = open(foz_to_end_blastout[f][j], "r")
            except:
                print "Could not open file: " + foz_to_end_blastout[f][j]
            
            lines = handle.readlines()
            handle.close()
            
            blast_hits = {}
            
            # convenience function to add hits
            def add_to_blast_hits(fields):
                if query not in blast_hits:
                    blast_hits[query] = {}
            
                dot_hits = dot_pattern.match(fields[1])
                if dot_hits:
                    lib = dot_hits.group(1)
                    forward_backward = dot_hits.group(2).split(".")[0]
            
                tail_hits = tail_pattern.match(fields[1])
                if tail_hits:
                    well = tail_hits.group(1)
                hit = lib + well
                if hit not in blast_hits[query]:
                    blast_hits[query][hit] = {}
                if forward_backward not in blast_hits[query][hit]:
                    blast_hits[query][hit][forward_backward] = []
                b_result = { "query": fields[0].strip(), 
                             "target": fields[1].strip(), 
                             "fb": forward_backward, 
                             "percent_sim": fields[2].strip(),
                             "a_length": fields[3].strip(), 
                             "q_start": int(fields[6].strip()),
                             "q_stop": int(fields[7].strip()),
                             "t_start": int(fields[8].strip()),
                             "t_stop": int(fields[9].strip()),
                             "score": int(fields[11].strip()),
                             "eval": fields[10].strip() }
                blast_hits[query][hit][forward_backward].append(b_result)
            
            for l in lines:
                fields = l.split("\t")
                query = fields[0].strip()
                q_start = int(fields[6].strip())
                q_end = int(fields[7].strip())
                # must hit near the full-lenth fosmid end
                length = int(xout_lengths[f][query]) - 10
                if q_start < q_end:
                    # forward
                    if q_start < 10 or q_end > length:
                        add_to_blast_hits(fields)
                else:
                    # backward
                    if q_end < 10 or q_start > length:
                        add_to_blast_hits(fields)
            
            if f not in final_mapping:
                final_mapping[f] = {}
            # collect all blast_hits
            for hit in blast_hits:
                if hit not in final_mapping[f]:
                    final_mapping[f][hit] = {}
                for k in blast_hits[hit]:
                    if k not in final_mapping[f][hit]:
                        final_mapping[f][hit][k] = {}
                    for l in blast_hits[hit][k]:
                        final_mapping[f][hit][k][l] = blast_hits[hit][k][l]
    return final_mapping
    
def rename_writeout_fosmids(rename_mapping, decomplex=True):
    # find targets in fasta
    
    org_fastas = get_files(temp_fasta_dir)
    
    for fasta in org_fastas:
        reader = FastaReader(temp_fasta_dir + sep + fasta)
        # open output_fasta
        try:
            fasta_out = fasta.rstrip(".xout")
            fasta_handle = open(output_dir + sep + fasta_out + ".qc", "w")
        except:
            "Problem opening fasta file"
        for record in reader:
            name = record.name.replace(">", "")
            seq = record.sequence
            if (fasta in rename_mapping) and (name in rename_mapping[fasta]):
                if len(rename_mapping[fasta][name]) > 0:
                    # it has a new name
                    new_name = rename_mapping[fasta][name]
                    name = "_".join(new_name)
                    # mels output of separating all found fosmids in their own file
                    if decomplex == True:
                        try:
                            t_fasta_out = fasta_out
                            mels_fasta = open(output_dir + sep + t_fasta_out + "." + name + ".fasta", "w")
                        except:
                            "Problem opening fasta file"
                        mels_fasta.write(">" + name + "\n")
                        mels_fasta.write(seq + "\n")
                
            fasta_handle.write(">" + name + "\n")
            fasta_handle.write(seq + "\n")
        fasta_handle.close()

def all_fosmidqc_report(final_mapping):
   results_handle = open(output_dir + sep + "all_fosmidqc_report" + ".txt", "w")
   line = ""
   line = "Fosmid QC Final Results:"
   print line
   results_handle.write(line + "\n")
   for f in xout_lengths:
       line = str(f) 
       f_temp = re.match("^(.*?)-",f)
       if f_temp:
           f_temp = f_temp.group(1)
           if f_temp in lib_name_maps['library_name_to_external']:
               line = line + " (" + str(lib_name_maps['library_name_to_external'][f_temp]) + ")"
       line = line + " (" + str(len(xout_lengths[f])) + " contig(s)):" 
       print line
       results_handle.write(line + "\n")
       for n in xout_lengths[f]:
           if len(xout_lengths[f]) > 1:
               # multiple fosmids print all hits
               line = ">" + str(n) + " ( " + str(xout_lengths[f][n]) + " nts ):"
               print line
               results_handle.write(line + "\n")
               if n not in final_mapping[f]:
                   line = "\tNo Hits"
                   print line
                   results_handle.write(line + "\n")
               else:
                   for e in final_mapping[f][n]:
                       line = "\t" + str(e)
                       print line
                       results_handle.write(line + "\n")
                       for r in final_mapping[f][n][e]:
                           line = "\t\t" + str(r)
                           print line
                           results_handle.write(line + "\n")
                           for j in final_mapping[f][n][e][r]:
                               line = "\t\t\t" + "(" + str(j["query"]) + "," + str(j["target"]) + "): " + " Length:" + str(j["a_length"]) + " Percent: " +  str(j["percent_sim"]) + " " + str(j["q_start"]) + "-" + str(j["q_stop"])
                               print line
                               results_handle.write(line + "\n")
           else:
               # only one fosmid

               line = ">" + str(n) + "( " + str(xout_lengths[f][n]) + " nts ):"
               print line
               results_handle.write(line + "\n")
               if n not in final_mapping[f]:
                   line = "\tNo Hits"
                   print line
                   results_handle.write(line + "\n")
               else:
                   for e in final_mapping[f][n]:
                       if len(final_mapping[f][n][e]) > 1:
                           line = "\t" + str(e)
                           print line
                           results_handle.write(line + "\n")
                           for r in final_mapping[f][n][e]:
                               line = "\t\t" + str(r)
                               print line
                               results_handle.write(line + "\n")
                               for j in final_mapping[f][n][e][r]:
                                   line = "\t\t\t" + "(" + str(j["query"]) + "," + str(j["target"]) + "): " + " Length:" + str(j["a_length"]) + " Percent:" + str(j["percent_sim"]) + " " + str(j["q_start"]) + "-" + str(j["q_stop"])
                                   print line
                                   results_handle.write(line + "\n")
                       elif len(final_mapping[f][n][e]) == 1:
                           line = "\t" + str(e)
                           print line
                           results_handle.write(line + "\n")
                           line = "\t**Only one hit**"
                           print line
                           results_handle.write(line + "\n")
                           for r in final_mapping[f][n][e]:
                               line = "\t\t" + str(r)
                               print line
                               results_handle.write(line + "\n")
                               for j in final_mapping[f][n][e][r]:
                                   line = "\t\t\t" + "(" + str(j["query"]) + "," + str(j["target"]) + "): " + " Length:" + str(j["a_length"]) + str(" Percent:") + str(j["percent_sim"]) + " " + str(j["q_start"]) + "-" +  str(j["q_stop"])
                                   print line
                                   results_handle.write(line + "\n")
   
   results_handle.close()
    


def main(argv): 
   args = vars(parser.parse_args())
   if check_arguments(args):
      print "ERROR: some required arguments not present or accessable"
      sys.exit(0)
      
   # set input parameters and global variables
   global input_dir
   global input_fastas
   global remove_me_dir
   global remove_me_fasta
   global blast_exe
   global blast_db_exe
   global naming_file_dir
   global naming_file
   global fosmid_end_dir
   global fosmid_ends
   global output_dir
   global xout_lengths
   
   input_dir = os.path.dirname(args['input_fastas'][0])
   input_fastas = paths_to_filenames(args['input_fastas'])
   remove_me_dir = os.path.dirname(args['remove_me_fasta'])
   remove_me_fasta = args['remove_me_fasta']
   blast_exe = args['blast_executable']
   blast_db_exe = args['database_executable']
   if args['naming_file']:
       naming_file_dir = os.path.dirname(args['naming_file'])
       naming_file = args['naming_file']
   fosmid_end_dir = os.path.dirname(args['fosmid_ends'][0])
   fosmid_ends = paths_to_filenames(args['fosmid_ends'])
   output_dir = args['output_directory']
   percent_id = args['percent_id']
   
   # create directories if they don't exist
   temp_dirs = [temp_blastdb_dir, temp_blastout_dir, temp_fasta_dir]
   for d in temp_dirs:
       if not os.path.exists(d):
           os.makedirs(d)
   if not os.path.exists(output_dir):
       os.makedirs(output_dir)
   
   # create lengths of fosmids and fosmid ends
   global fosmid_lengths
   for f in input_fastas:
       if f not in fosmid_lengths:
           fosmid_lengths[f] = {}
       fosmid_lengths[f] = fasta_lengths(input_dir + sep + f)   
   
   # parse the global lib name mapping file
   if args['naming_file']:       
       global lib_name_maps
       lib_name_maps = parse_naming_file(naming_file)
   
   # create the remove_me_database
   xout_db_location = create_database(remove_me_fasta, blast_db_exe, temp_blastdb_dir)
   
   # for each input fasta file, blast against the target database
   fasta_to_blastout_remove = {}
   for f in input_fastas:
       blastout_file = blast(blast_exe, input_dir + sep + f, xout_db_location, percent_id, temp_blastout_dir)
       # pair fasta with blastout file
       fasta_to_blastout_remove[f] = blastout_file
   
   # finish the xout step
   remove_primer_sequences(fasta_to_blastout_remove, int(args['min_length']))
   
   if args['just_qc']:
       exit()
   
   # clean up temp_blastout, temp_blastdb
   clean_dir(temp_blastdb_dir)
   clean_dir(temp_blastout_dir)
   
   # get all xout files
   global xout_lengths
   xout_files = get_files(temp_fasta_dir)
   
   for x_file in xout_files:
       if x_file not in xout_lengths:
           xout_lengths[x_file] = {}
       xout_lengths[x_file] = fasta_lengths(temp_fasta_dir + sep + x_file)
   
   # create blast database for fosmid ends
   f_end_to_db_loc = {} # fos_end db lookup
   for f_end in fosmid_ends:
       db_location = create_database(fosmid_end_dir + sep + f_end, blast_db_exe, temp_blastdb_dir)
       f_end_to_db_loc[f_end] = db_location
   
   # qc_fosmid to end
   foz_to_end_blastout = {}
   for qc_file in xout_files:
       if qc_file not in foz_to_end_blastout:
           foz_to_end_blastout[qc_file] = {}
       for f_end in fosmid_ends:
           if f_end not in foz_to_end_blastout[qc_file]:
               foz_to_end_blastout[qc_file][f_end] = {}
           blastout_file = blast(blast_exe, temp_fasta_dir + sep + qc_file, f_end_to_db_loc[f_end], percent_id, temp_blastout_dir)
           foz_to_end_blastout[qc_file][f_end] = blastout_file
   
   
   final_mapping = identify_fosmid_matches(foz_to_end_blastout)
   
   all_fosmidqc_report(final_mapping)
   
   rename_mapping = {}
   
   fb_handle = open(output_dir + sep + "fosmidqc_mapping_report" + ".txt", "w")
   
   line = "Fosmid QC Forward-Backward Hit Mapping:"
   print line
   fb_handle.write(line + "\n")
   count = 0
   for i in xout_lengths:
       line = i
       f_temp = re.match("^(.*?)-",i)
       if f_temp:
          f_temp = f_temp.group(1)
          if f_temp in lib_name_maps['library_name_to_external']:
              line = line + " (" + str(lib_name_maps['library_name_to_external'][f_temp]) + ")"
       print line
       fb_handle.write(line + "\n")
       if i not in rename_mapping:
           rename_mapping[i] = {}
       for j in xout_lengths[i]:
           if j not in rename_mapping[i]:
               rename_mapping[i][j] = []
           if j in final_mapping[i]:
               for k in final_mapping[i][j]:
                   if len(final_mapping[i][j][k]) > 1:
                       # found both ends append to list
                       rename_mapping[i][j].append(k)
           # print mapping
           if len(rename_mapping[i][j]) > 0:
               line = "\t" + ">" + str(j) + " --> " + str(rename_mapping[i][j])
               print line
               fb_handle.write(line + "\n")
               count += 1
           else:
               line = "\t" + ">" + str(j) + " --> " + "None"
               print line
   
   fb_handle.close()
   # create qc_fosmids
   rename_writeout_fosmids(rename_mapping)
   manual_aid = open(output_dir + sep + "fosmidqc_manual_aid_report" + ".txt", "w")
   line = "Fosmid QC Manual Aid Report:"
   rename_mapping = {}
   manual_aid.write(line + "\n")
   for i in xout_lengths:
       if i not in rename_mapping:
           rename_mapping[i] = {}
       for j in xout_lengths[i]:
          if j not in rename_mapping[i]:
              rename_mapping[i][j] = []
          if j in final_mapping[i]:
              for k in final_mapping[i][j]:
                  if len(final_mapping[i][j][k]) == 1:
                      # found both ends append to list
                      rename_mapping[i][j].append(k)
          # # print mapping
          # if len(rename_mapping[i][j]) > 1:
          #     line = "\t" + ">" + str(j) + " --> " + str(rename_mapping[i][j])
   
   for i in rename_mapping:
       for j in rename_mapping[i]:
           if len(rename_mapping[i][j]) > 0:
               mapping = ""
               f_temp = re.match("^(.*?)-",i)
               if f_temp:
                  f_temp = f_temp.group(1)
                  if f_temp in lib_name_maps['library_name_to_external']:
                     mapping = " (" + str(lib_name_maps['library_name_to_external'][f_temp]) + ") "
               line = str(i) + mapping + ">" + str(j) + " " + str(rename_mapping[i][j])
               manual_aid.write(line + "\n")
   manual_aid.close()
   
   print str(count) + " Fosmid sequence(s) Forward-Backward Hit mapped."
   
   # remove temporary directories
   for i in temp_dirs:
      clean_dir(i)
      os.removedirs(i)
   
   exit()
   


# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
