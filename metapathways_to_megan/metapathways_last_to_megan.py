#!/usr/bin/python

import sys
import re

# usage: python metapathways_last_to_megan.py -i input_filename -o output_filename

input_filename = sys.argv[2]
output_filename = sys.argv[4]
file_handle = open(input_filename, "r")

lines = file_handle.readlines()

file_handle.close()

brackets_pattern = re.compile("\[(.*)\]")

output_handle = open(output_filename, "w")

for l in lines:
	fields = l.split("\t")
	hits = brackets_pattern.search(fields[8])
	if hits:
	   read = fields[0]
	   last_score = fields[2]
	   taxa = hits.group(1)
	   out_line = read + ", " + taxa + ", " + last_score + "\n"
	   output_handle.write(out_line)

output_handle.close()

exit()