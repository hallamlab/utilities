#!usr/bin/python
# Takes correlation table and generates network files (nodes & edges).


try:
     import os
     import re
     import sys
     import argparse
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)


# Example: python metapathways_last_to_megan.py -i output/HOT_Sanger/*/blast_results/*cog*out.txt -o .
what_i_do = """Retrives read-taxonomy hits from the parsed (blast/last) result files <sample>.<database>.(blast/last)out.parsed.txt files, 
e.g. my_sample.refseq.lastout.parsed.txt, files from the <sample>/blast_results/ directory and formats them as .csv 
files for import into MEGAN. Requires that the database keep the taxonomy in within square braces, e.g. [E. coli K12].
"""
parser = argparse.ArgumentParser(description=what_i_do)
# add arguments to the parser
parser.add_argument('-i', dest='input_cor', type=str, nargs='+',
                required=True, help='input correlation matrix is csv format', default=None)          
parser.add_argument('-o', dest='output_dir', type=str, nargs='?',
                required=False, help='result directory where results are placed', default=os.getcwd())
parser.add_argument('--full', dest='full', action='store_true',
                required=False, help='flag to calculate the full matrix', default=False)              
parser.add_argument('--cor_cutoff', dest='cor_cutoff', type=str, nargs='?',
                required=False, help='absolute correlation cuttoff (default |0.3|)', default=0.3)
parser.add_argument('--pval', dest='pval_file', type=str, nargs='?',
                    required=False, help='corresponding matrix of pvalues', default=None)
parser.add_argument('--pval_cutoff', dest='pval_cutoff', type=str, nargs='?',
                    required=False, help='pval cuttoff (default <0.05)', default=0.05)
                
def main(argv):
    # get arguments
    args = vars(parser.parse_args())
    
    input_cor = args["input_cor"]
    out_dir = args["output_dir"]
    
    # script could take many files
    for f in input_cor:
        file_handle = open(f, "r")
        lines = file_handle.readlines()
        header = lines[0].split("\t")
        
        edges_out = open(os.path.abspath(args["output_dir"]) + os.sep + os.path.basename(f) + "_edges.txt", "w")
        nodes_out = open(os.path.abspath(args["output_dir"]) + os.sep + os.path.basename(f) + "_nodes.txt", "w")
        
        nodes_out.write("node_name" + "\n")
        # chris wants all the nodes
        for n in header[1:len(header)]:
            node_line = n + "\n"
            nodes_out.write(node_line)
        
        nodes_out.close() # close node file
        
        # write the edge header
        edge_header = "from" + "\t" + "to" + "\t" + "corr"
        if args["pval_file"]:
           edge_header = edge_header + "\t" + "pval"
        
        edge_header = edge_header + "\n"
        edges_out.write(edge_header)
        
        pval_lines = None
        
        if args["pval_file"]:
            pval_file_handle = open(args["pval_file"], "r")
            pval_lines = pval_file_handle.readlines()
            pval_file_handle.close()
        
        # for each line in the file
        for i in range(1,len(lines)):
            fields = lines[i].split("\t")
            if args["pval_file"]:
                p_fields = pval_lines[i].split("\t")
            if args["full"]:
                # creates an edge for all values in the matrix
                for j in range(1,len((fields))):
                    if abs(float(fields[j].strip("\n").strip())) > abs(float(args["cor_cutoff"] )):
                        edge_out_line = fields[0].strip("\n").strip() + "\t" + \
                                        header[j].strip("\n").strip() + "\t" + \
                                        fields[j].strip("\n").strip()
                        if args["pval_file"]:
                            if (float(p_fields[j].strip("\n").strip()) <= float(args["pval_cutoff"])):
                                edge_out_line = edge_out_line + "\t" + p_fields[j].strip("\n").strip()
                            else:
                                # bail!
                                continue
                        edge_out_line = edge_out_line + "\n"
                        edges_out.write(edge_out_line)
            else:
                # just calculate the upper triangular
                for j in range(i, len((fields))):
                    if abs(float(fields[j].strip("\n").strip())) > abs(float(args["cor_cutoff"] )):
                        edge_out_line = fields[0].strip("\n").strip() + "\t" + \
                            header[j].strip("\n").strip() + "\t" + \
                            fields[j].strip("\n").strip()
                        if args["pval_file"]:
                            if (float(p_fields[j].strip("\n").strip()) <= float(args["pval_cutoff"])):
                                edge_out_line = edge_out_line + "\t" + p_fields[j].strip("\n").strip()
                            else:
                                # bail!
                                continue
                        edge_out_line = edge_out_line + "\n"
                        edges_out.write(edge_out_line)

        # close edges file
        edges_out.close()
    
    exit()
    
# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])







