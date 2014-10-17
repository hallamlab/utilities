assembly_stats.py
-----------------

A simple utility to calculate assembly statistics from `.fasta` files.

## Example Usage

```
python assembly_stats.py -i <fasta GLOB> -o <output_file.txt> [--minl] <min_length>
```
where,

* `-i`: glob to select fasta files (e.g., `*.fasta`)
* `-o`: name of output file
* `--minl`: minimum length cutoff

### Example Output

* script outputs a tab-delimited file

```
file	stat	value
Sak_2011_05_24_120m_454_Contigs_500.fasta	N	30217
Sak_2011_05_24_120m_454_Contigs_500.fasta	N_Trimmed	30217
Sak_2011_05_24_120m_454_Contigs_500.fasta	Sum_l	24971863
Sak_2011_05_24_120m_454_Contigs_500.fasta	Min_l	501
Sak_2011_05_24_120m_454_Contigs_500.fasta	Max_l	21944
Sak_2011_05_24_120m_454_Contigs_500.fasta	Med_l	680
Sak_2011_05_24_120m_454_Contigs_500.fasta	Mean_l	826
Sak_2011_05_24_120m_454_Contigs_500.fasta	N5	3719
Sak_2011_05_24_120m_454_Contigs_500.fasta	N10	2439
Sak_2011_05_24_120m_454_Contigs_500.fasta	N15	1796
Sak_2011_05_24_120m_454_Contigs_500.fasta	N20	1393
Sak_2011_05_24_120m_454_Contigs_500.fasta	N25	1164
```
