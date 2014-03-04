#!/usr/bin/python

from __future__ import division
try:
    import re
    import sys
    import os
    import types
    import argparse
    
except:
    import_exception = """ Could not load some modules """
    print import_exception;
    
what_i_do = "inserts a list of organisms into the NCBI Newick format"

parser = argparse.ArgumentParser(description=what_i_do)
parser.add_argument('--tree', dest='tree', type=str, nargs=1,
                required=True, help='the input tree in newick format', default=None)
parser.add_argument('--id_list', dest='id_list', type=str, nargs=1,
                required=False, help='list of translation ids for NCBI names', default=None)                
parser.add_argument('--list', dest='list', type=str, nargs=1,
               required=True, help='a list of parent child relations by id', default=None)
parser.add_argument('--out', dest='out', type=str, nargs=1,
               required=True, help='output file name of updated tree', default=None)

def fnd(fname, s, start=0):
    with open(fname, 'rb') as f:
        fsize = os.path.getsize(fname)
        bsize = 1
        buffer = None
        if start > 0:
            f.seek(start)
        overlap = len(s) - 1
        while True:
            if (f.tell() >= overlap and f.tell() < fsize):
                f.seek(f.tell() - overlap)
            buffer = f.read(bsize)
            if buffer:
                pos = buffer.find(s)
                if pos >= 0:
                    return f.tell() - (len(buffer) - pos)
            else:
                return -1

  
# the main function
def main():
  # parse arguments
  args = vars(parser.parse_args())
  
  id_to_name = {}
  name_to_id = {}
  
  # print "Loading ID Map:"
  # id_handle = open(args['id_list'][0], 'r')
  # id_lines = id_handle.readlines()
  # id_handle.close()
  # 
  # for l in id_lines:
  #     fields = l.split("\t")
  #     name_to_id [ fields[1].strip() ] = fields[0].strip()
  #     id_to_name [ fields[0].strip() ] = fields[1].strip()
  # print "Done."
  
  print "Loading New Children List:"
  list_handle = open(args['list'][0], 'r')
  parent_to_new = {} # parent to list of children to add
  positions = [] # list of parents to add children to
  
  # create the list of children to add
  for l in list_handle:
      fields = l.split("\t")
      parent = fields[0].strip()
      child = fields[1].strip()
      if parent not in parent_to_new:
          # create child list
          parent_to_new[parent] = []
          positions.append(parent)
      parent_to_new[parent].append(child)
  list_handle.close()
  print "Done."
  
  print "Opening Tree:"
  tree_buff = open(args['tree'][0], "rb")
  
  num_flag = False
  num = ""
  char_old = None
  char = None
  tree_buff.seek(0, os.SEEK_END)
  length = tree_buff.tell()
  print "  - scanning and adding new children to tree"
  temp_out = open("rev", "w")
  for i in range(length, 0, -1):
      tree_buff.seek(i-1,0)
      char = tree_buff.read(1)
      if (re.match(r'[0-9]', char)):
          # is a number
          if not num_flag:
             # reset number
             num = ""
             num_flag = True
          num = char + num
      else:
          if num_flag:
              # string[::-1] # python way to reverse a string
              # print num
              # check if number equals our search candiate
              if num in positions:
                  if char == "(" or char == ",":
                      # no children: print to the left with parenthesis
                      l = char + "(" + ",".join(parent_to_new[num]) + ")" + num
                      # print l
                      temp_out.write(l[::-1])
                  elif char == ")":
                      # children: place into children
                      l = ""
                      for c in parent_to_new[num]:
                          l = l + "," + str(c)
                      l = l + char + num
                      # print l
                      temp_out.write(l[::-1])
                  num = ""
                  num_flag = False
              else:
                  l = char + num
                  temp_out.write(l[::-1])
                  num = ""
                  num_flag = False
          else:
              # reset the number
              temp_out.write(char)
              num = ""
              num_flag = False
  
  # TODO: test this corner case
  # if num in positions:
  #     print num
  
  tree_buff.close()
  temp_out.close()
  print "Done."
  
  # write out reverse tree to disk as rev and then
  # reverse into output file for correct order.
  print "Writing out new tree."
  rev_handle = open("rev", "rb")
  out_handle = open(args['out'][0], "w")
  
  rev_handle.seek(0, os.SEEK_END) # get position of end of file
  length = rev_handle.tell()
  for i in range(length, 0, -1):
      rev_handle.seek(i-1,0)
      char = rev_handle.read(1)
      out_handle.write(char)
  rev_handle.close()
  out_handle.close()
  
  # remove the reversed file
  os.remove("rev")
  print "Done."
  
  exit()

# call the main function
if __name__ == "__main__":
   sys.exit(main())