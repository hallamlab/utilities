#!/usr/bin/python

try:
     import os
     import re
     import sys
     import argparse
except:
    print "Could not load some modules"

what_i_do = "Given a list of fasta files, calculates a series of assembly statistics in longtable format"
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='inputfastas', type=str, nargs='+',
                required=True, help='a selection of one or many fosmid input .fasta files (required)', default=None)
parser.add_argument('-o', dest='output', type=str, nargs='?',
                required=False, help='output file (required)', default="assembly_stats.txt")
parser.add_argument('--minl', dest='minlength', type=int, nargs='?',
                required=False, help='a selection of one or many fosmid input .fasta files (required)', default=50)

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


def nx(lens, percent):
    """
       Calculates any NXX (e.g. N50, N90) statistic.
    """
    
    lenSum = sum(lens)
    threshold = (float(percent) / 100) * lenSum
    
    runningSum = 0
    nxx = 0
    nxxLen = 0

    for i in range(len(lens)-1, -1, -1):
       myLen = lens[i]
       nxx += 1
       runningSum += myLen

       if runningSum >= threshold:
          nxxLen = myLen
          break

    return nxx, nxxLen

def trimLens(lens, minLen):
   """
   Eliminates any reads below a certain threshold.  Function assumes that input
   list lens is sorted smallest to largest.
   """
   
   index = 0
   
   for i in range(len(lens)):
      if lens[i] < minLen:
         index += 1
      else:
         break
   
   return lens[index:len(lens)]


def getLens(filename):
   """
   Parses FASTA file using screed to create a sorted list of contig lengths.
   """
   lens = []
   reader = FastaReader(filename)
   
   for record in reader:
       lens.append(len(record.sequence))
   
   return sorted(lens)


def main(argv):
    """
    Main function.
    """
    
    # parse arguments
    args = vars(parser.parse_args())
    
    # list of fastas
    fastas = args['inputfastas']
    minlength = args['minlength']
    
    # output
    stats_out = open(args['output'],"w")
    
    header = "file\tstat\tvalue\n"
    stats_out.write(header)
    
    # for each filename calculate some stats
    for filename in fastas:
        lens = getLens(filename) # list of lengths
        trimmedLens = trimLens(lens, minlength) #
        
        if len(trimmedLens) == 0:
            print filename + " - no sequences longer than threshold"
            continue
        
        statN = len(lens)
        statTrimmedN = len(trimmedLens)
        statSum = sum(trimmedLens)
        statMin = min(trimmedLens)
        statMax = max(trimmedLens)
        statMed = trimmedLens[ (len(trimmedLens)-1) / 2 ]
        statMean = int( statSum / float(statTrimmedN) )
        statnxs = []
        statnx_lens = []
        for n in range(5,96,5):
            statnx, statnx_len = nx(trimmedLens, n)
            # print statnx, statnx_len
            statnxs.append("N" + str(n))
            statnx_lens.append(statnx_len)
        
        stats = [statN, statTrimmedN, statSum, statMin, statMax, statMed, statMean] + list(statnx_lens)
        stat_names = ["N", "N_Trimmed", "Sum_l", "Min_l", "Max_l", "Med_l", "Mean_l"] + statnxs
        
        for i in range(len(stats)):
            lin= filename + "\t" + str(stat_names[i]) + "\t" + str(stats[i])  + "\n"
            stats_out.write(lin)
        
        # also writeout lengths of each file
        lengths_handle = open(filename + ".lens.txt", "w")
        trimmedLens.append("\n")
        lengths_handle.write("\t".join(map(str, trimmedLens)))
        lengths_handle.close()
    
    stats_out.close()
    
    



# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])