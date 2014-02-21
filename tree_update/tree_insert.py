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
                required=False, help='the input tree in newick format', default=None)
parser.add_argument('--id_list', dest='id_list', type=str, nargs=1,
                required=False, help='list of translation ids for NCBI names', default=None)                
parser.add_argument('--list', dest='list', type=str, nargs=1,
               required=False, help='a list of parent child relations by id', default=None)
parser.add_argument('--out', dest='out', type=str, nargs=1,
               required=False, help='output names', default=None)

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
  
  print "Loading ID Map:"
  id_handle = open(args['id_list'][0], 'r')
  id_lines = id_handle.readlines()
  id_handle.close()
  
  for l in id_lines:
      fields = l.split("\t")
      name_to_id [ fields[1].strip() ] = fields[0].strip()
      id_to_name [ fields[0].strip() ] = fields[1].strip()
  print "Done."
  
  
  print "Loading List:"
  list_handle = open(args['list'][0], 'r')
  parent_to_new = {}
  for l in list_handle:
      fields = l.split("\t")
      parent_to_new[fields[0].strip()] = fields[1].strip()
  
  list_handle.close()
  print "Done."
  positions = []
  for key, value in parent_to_new.iteritems() :
     positions.append(key)
  
  print "Opening Tree:"
  tree_buff = open(args['tree'][0], "rb")
  
  num_flag = False
  num = ""
  char_old = None
  char = None
  for i in range(length, 0, -1):
      tree_buff.seek(i-1,0)
      char = tree_buff.read(1)
      if (re.match(r'[0-9]', char)):
          # is a number
          if not num_flag:
             # reset number
             num = ""
             num_flag = True
          num = num + char
      else:
          if new_num:
              # check if number equals our search candiate
              if num in positions:
                  if char_old == "(" or char_old == ",":
                      # no children: print to the left with parenthesis
                      l = "(" + str(parent_to_new[num]) + ")" + num
                      sys.stdout.write(l)
                  elif char_old == ")":
                      # children: place into children
                      l = "," + str(parent_to_new[num]) + ")" + num
                      sys.stdout.write(l)
                  # print tree_buff.tell()
              else:
                  l = num
                  sys.stdout.write(l)
          # reset the number
          background = background + char
          num = ""
          num_flag = False
          char_old = char # remember 
  if num in positions:
      print num
      print tree_buff.tell()
  print "Done."
  
  # for p in positions:
  #     tree_buff.seek(p-10)
  #     print tree_buff.read(20)
  
  exit()

# call the main function
if __name__ == "__main__":
   sys.exit(main())