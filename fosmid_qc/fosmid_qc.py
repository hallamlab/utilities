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
     from os import makedirs, sys, remove
     from sys import path
     from libs.metapaths_utils  import parse_command_line_parameters, fprintf
     from libs.sysutil import getstatusoutput, pathDelim
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)

what_i_do = "I'm a script that does something"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_fastas', type=str, nargs='+',
                required=True, help='a selection of one or many input .fasta files (required)', default=None)
parser.add_argument('-d', dest='remove_me_fasta', type=str, nargs='?',
                required=True, help='a fasta file containing sequences to find and remove from input fastas (required)', default=None)                
parser.add_argument('-o', dest='output_directory', type=str, nargs='?',
                required=True, help='the name a target output directory (required)', default=None)
parser.add_argument('-b', dest='blast_executable', type=str, nargs='?',
                required=False, help='location of the blastn executable, will assume in PATH if not specified', default='blastn')
parser.add_argument('-f', dest='database_executable', type=str, nargs='?',
                required=False, help='location of the makeblastdb executable, will assume in PATH if not specified', default='makeblastdb')
               
def check_arguments(args):
   if args['input_fastas'] == None or args['output_directory'] == None or args['remove_me_fasta'] == None:
      # TODO check to see if these files/locations are actually accessable
      return True
   else:
      return False

class FastaRecord(object):
   def __init__(self, name, sequence):
       self.name = name
       self.sequence = sequence

#    return FastaRecord(title, sequence)

def read_fasta_records(input_file):
   records = []
   sequence=""
   name=""
   while 1:
        line = input_file.readline()
        if line == "": 
           if sequence!="" and name!="":
              records.append(FastaRecord(name, sequence))
           return  records

        if line=='\n':
           continue

        line = line.rstrip()
        if  line.startswith(">") :
           if sequence!="" and name!="":
              records.append(FastaRecord(name, sequence))

           name = line.rstrip()
           sequence =""
        else:
           sequence = sequence + line.rstrip()
   return records

def format_db_blast(formatdb_executable, seq_subset_file):
   cmd='%s -dbtype prot -in %s' %(formatdb_executable, seq_subset_file.name)
   result= getstatusoutput(cmd)


def format_db_last(formatdb_executable, seq_subset_file):
   dirname = os.path.dirname(seq_subset_file.name)     
   cmd='%s -p -c %s  %s' %(formatdb_executable, dirname + PATHDELIM + 'subset_db', seq_subset_file.name)
   result= getstatusoutput(cmd)


def blast_against_itself(blast_executable, seq_subset_file, blast_table_out):
   cmd='%s -outfmt 6 -db  %s -query %s -out  %s' %(blast_executable,  seq_subset_file.name, seq_subset_file.name, blast_table_out)
   result= getstatusoutput(cmd)

def last_against_itself(last_executable, seq_subset_file, last_table_out):
   dirname = os.path.dirname(seq_subset_file.name)     
   cmd='%s -o %s -f 0 %s %s' %(last_executable,  last_table_out, dirname + PATHDELIM + 'subset_db',  seq_subset_file.name)
   result= getstatusoutput(cmd)


def add_last_refscore_to_file(blast_table_out, refscore_file, allNames):
   commentPATTERN = re.compile(r'^#')

   infile = open( blast_table_out,'r')
   refscores = {}
   lines = infile.readlines()
   for line in lines:
      if commentPATTERN.match(line):
         continue
      line=line.rstrip()
      fields = line.split('\t')
      if len(fields) != 12:
         print 'Error in the blastout file'
         sys.exit(1)
      if fields[6].rstrip()==fields[1].rstrip():
     #    fprintf(refscore_file, "%s\t%s\n",fields[0], fields[11])
         refscores[fields[1]]=fields[0]

   for key, value in refscores.iteritems():
      allNames[key] = True
      fprintf(refscore_file, "%s\t%s\n",key, value)

   infile.close()



def add_blast_refscore_to_file(blast_table_out, refscore_file, allNames):
   infile = open( blast_table_out,'r')
   refscores = {}
   lines = infile.readlines()
   for line in lines:
      line=line.rstrip()
      fields = line.split('\t')
      if len(fields) != 12:
         print 'Error in the blastout file'
         sys.exit(1)
      if fields[0].rstrip()==fields[1].rstrip():
     #    fprintf(refscore_file, "%s\t%s\n",fields[0], fields[11])
         refscores[fields[0]]=fields[11]

   for key, value in refscores.iteritems():
      allNames[key] = True
      fprintf(refscore_file, "%s\t%s\n",key, value)

   infile.close()


# compute the refscores
def compute_refscores(formatdb_executable, blast_executable,seq_subset_file, refscore_file, allNames, algorithm):
   if algorithm =='LAST':
       format_db_last(formatdb_executable, seq_subset_file)
       last_table_out = seq_subset_file.name + ".lastout"
       last_against_itself(blast_executable, seq_subset_file, last_table_out)
       add_last_refscore_to_file(last_table_out,refscore_file, allNames)

   if algorithm =='BLAST':
      format_db_blast(formatdb_executable, seq_subset_file)
      blast_table_out = seq_subset_file.name + ".blastout"
      blast_against_itself(blast_executable, seq_subset_file, blast_table_out)
      add_blast_refscore_to_file(blast_table_out,refscore_file, allNames)
   
   return None

def remove_blast_index_files(filename):
   prefixes = [ 'blastout', 'phr', 'pin', 'psq' ] 
   for prefix in prefixes:
      try:
         remove(filename +"." + prefix)
      except IOError:
         pass


def remove_last_index_files(filename):
   prefixes = [ 'des', 'sds', 'suf', 'bck',  'ssp', 'tis' ]

   remove( filename+ '.lastout')
   dirname = os.path.dirname(filename)     
   remove( dirname + PATHDELIM + 'subset_db' +'.prj')
   for prefix in prefixes:
      try:
         remove(dirname + PATHDELIM + 'subset_db0.' + prefix)
      except IOError:
         pass

def create_database(in_fasta, exe):
    print "in create_database"
    cmd='%s -dbtype nucl -in %s' %(exe, in_fasta)
    print cmd
    exit()
    result= getstatusoutput(cmd)
   

# the main function
SIZE = 1000

def main(argv): 
   args = vars(parser.parse_args())
   if check_arguments(args):
      print "ERROR: some required arguments not present or accessable"
      sys.exit(0)

   input_fastas = args['input_fastas']
   output_directory = args['output_directory']
   remove_me_fasta = args['remove_me_fasta']
   blast_executable = args['blast_executable']
   database_executable = args['database_executable']
   
   # create the absolute paths of all inputfiles
   for i in range(len(input_fastas)):
       input_fastas[i] = os.path.abspath(input_fastas[i])
   output_directory = os.path.abspath(output_directory)
   remove_me_fasta = os.path.abspath(remove_me_fasta)
   blast_executable = blast_executable
   database_executable = database_executable
   
   print input_fastas
   print output_directory
   print remove_me_fasta
   print blast_executable
   print database_executable
   
   # create the remove_me_database
   create_database(remove_me_fasta, database_executable)
   exit()

   # input file to blast with itself to commpute refscore
   infile = open(input_fasta,'r')

   #this file has the refscores of the entire file
   outfile = open(output_fasta, 'w') 

   count = 0

   allNames= dict()
   for record in read_fasta_records(infile):
       if count % SIZE == 0:
           if count > 0:
             seq_subset_file.close()
             compute_refscores(formatdb_executable, blast_executable,seq_subset_file, outfile, allNames, algorithm);

             # now remove the old file
             if algorithm == 'BLAST' :
                remove_blast_index_files(seq_subset_file.name)

             if algorithm == 'LAST' :
                remove_last_index_files(seq_subset_file.name)

             remove(seq_subset_file.name)

           seq_subset_file = open(output_file +'.tmp.'+ str(count) +'.fasta','w')
       allNames[record.name.replace(">","")] = False;    
       fprintf(seq_subset_file, "%s\n", record.name)
       fprintf(seq_subset_file, "%s\n", record.sequence)

       count = count + 1

   #print str(count) + "   "  + "going to blast last sequence "
   if (count) % SIZE != 0:
      #print str(count) + "   "  + "last sequence "
      seq_subset_file.close()
      compute_refscores(formatdb_executable, blast_executable,seq_subset_file, outfile, allNames, algorithm);
      remove(seq_subset_file.name)
      if algorithm == 'BLAST' :
         remove_blast_index_files(seq_subset_file.name)
      if algorithm == 'LAST' :
         remove_last_index_files(seq_subset_file.name)


   #print count
   for key in allNames:
       if allNames[key] ==False:
          fprintf(outfile, "%s\t%s\n",key, 1000000)
   
   outfile.close()

# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])

