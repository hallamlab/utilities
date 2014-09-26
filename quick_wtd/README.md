quick_wtd
---------

A quick utility script for calculating WTD on legacy output of MetaPathways.

This will produce a pathway table simmilar to `extract_pathways_from_pgdb.pl` included
with [MetaPathways v1.0](https://github.com/hallamlab/MetaPathways/), but addes columns
for the Weighted Taxonomic Distance (WTD), the Observed Pathway LCA Taxonomy, and closest member
of the particular pathways Expected Taxonomic Range(s).

The script requires:

* a `functional_and_taxonomic_table.txt` mapping annotated ORFs to Lowest Common Ancestor
 (LCA) taxonomies
* an ePGDB to be made with the ORFs named in the `functional_and_taxonomic_table.txt` file
* the `ncbi.map` file of MEGAN preferred taxonomic names for easier interpretation
* the `ncbi_taxonomy_tree.txt` containg a list of nodes in the NCBI Taxonomy Database

## Basic Useage:

`MetaPathways_extract_pathways.py` accepts the following options:

* `-o`: output filename
* `-p`: name of the pgdb
* `-t`: location of the pathway tools executable
* `-wtd`: flag to indicate the calculation of the WTD
* `--annotation-table`: location of the `functional_and_taxonomic_table.txt` file
* `--ncbi-tree`: location fo the `ncbi_taxonomy_tree.txt` file
* `--ncbi-megan-map`: location of the `ncbi.map` file

#### Example:
```
python MetaPathways_extract_pathways.py -o /tmp/test_table.txt \
                                        -p e4093112_combined_unique \
                                        -t /Users/nielsh/pathway-tools-17/pathway-tools/pathway-tools \
                                        --wtd \
                                        --annotation-table resources/functional_and_taxonomic_table.txt \
                                        --ncbi-tree resources/ncbi_taxonomy_tree.txt \
                                        --ncbi-megan-map resources/ncbi.map
```
