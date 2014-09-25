**********************************************************************
         CORRELATION NETWORK CALCULATOR


The "corrlation_network.pl" is a perl script that can be used to compute correlation
network from a OTU table (such as the one produced by QIIME). The requirements for 
running this script is the availability of a perl interpreter and the perl package 
Statistics::Distributeions (http://search.cpan.org/~mikek/Statistics-Distributions-1.02/Distributions.pm)


The script produces two files (i) the node file, which contains the list of 
nodes and their attributes, and the (ii) the edge file, which lists all the edges
in the network with the cutoffs provided. 

The script can be used as
Usage: ./correlation_network.pl -i inputfile   -o outputfile  

OPTIONS:
               [-min_count c (default 25%) ]  [-mlabel label ]  [ -max_negative_correl 0.x (default -0.4) ] 
               [ -min_positive_correl 0.y (default 0.4) ] [-max_p_value max_pvalue (default 0.05] 
               [ -min_threshold x] 

               This program will prepare the properties of the node 
                  --inputfile the otu table file
                  --outputfile is the name of the output file which produces two files: an edge file and a node file 
                  --labelfile  a label file that has the name of the sub samples, one in each row
                  --min_count  min number (in %)  of non-zero columns for an otu to be used, default is 25% of the total columns
                  --min_threshold  is the minimum normalized value of the read counts to count as non-zero

                 This script will create the edge and node files for the network
                  --statistic  [Pearson | Kendall | Spearman ] (Default Pearson)
                  --shuffle [ shuffles the columns of values row-by-row (default 10)  ] 


"inputfile" is the name of the OTU table file with tab delimited columns. The "outputfile" is 
a prefix appended to the node file and the edge file. For example, if "outputfile" is 
 "abc" then  the node and the edge files are named as "abc_nodes.txt" and "abc_edges.txt", respectively.

The corrleations can be computed for Person, Kendall and Sperman rank correlation with suppfiled p-value 
cutoff.  In addition the options "-min_positive_correlation" and "-max_negative_correl" are the 
limits of the correlation values that are deemed as significant correlation for the purpose of 
the network. 

Once the node and the edge files are created they can be imported to Cytoscape to visualize the network.
Severla different network visulation layouts are available in Cytoscape, one of the interesting and
relevant ones is the "force directed layout".
