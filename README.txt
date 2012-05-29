David L. Stern
Janelia Farm Research Campus
29 May 1012

This program finds oligonucleotide probe sequencing suitable for 
detecting polymorphisms using gDNA hybridization to probes on glass slides.
For each position, one probe is made to straddle the deletion and one is 
specific to the insertion.
Simple screens are performed to exclude runs of NNs and s1ple repeat regions.

DEPENDENCIES

Python 2.7
BioPython

USAGE

(If necessary, make file executable)
chmod -x indel_probes.1.5

./indel_probes.1.5 


Two items required for input:
(1) a folder of files, each containing a base-for-base sequence 
alignment of genomic regions from two strains/species in fasta format.

e.g.

>10.0_10.1_e A new nucleotide sequence entered manually
CGAAGCACCTGACGGATTTCTTGGCTGGTCCATGATTAGATAAATGGAAGCCTTTCT
>10.0_10.1_i A new nucleotide sequence entered manually
CGAAGCACCTGACTGATTTCTTGGCTGGTCCATGATTAGATAAATGGAAGCCTTTCT


(2) a 'map.txt' file of four tab-delimited columns that provides information 
about the alignments. (Note: file_name must match file names precisely, chromosome must be 'Chr..')

file_name chromosome start_position stop_position

e.g.

1	ChrX	10000000	10100000
2	ChrX	10100001	10200000
3	ChrX	10200005	10300000
3	ChrX	10300002	10400000

The output will be a list of probes for the two strains/species with each probe
indicated by 

alignmentName_chromosome_postion_<in or del>

Note, only the first contiguous string in the alignment name is used, so this should include a 
unique strain/species identified

e.g. 

10.0_10.1_e_ChrX_100003401_d	CCTAACTCCAACTACAACGACGATTTGATGAGTGCAGAGTACATAACAAACTAACGCCCC
10.0_10.1_i_ChrX_100003401_i	TACACCTAACTCCAACTACAACGACGAGACACTTAGGAATGTTTGATGAGTGCAGAGTCC
10.0_10.1_e_ChrX_100024401_d	TGGTCAGAAAACCAAGATCATCAACCCGTGTTATTCTGTCCAAAGGGCACCTACTTCGAG
10.0_10.1_i_ChrX_100024401_i	CGAAATATGGTCAGAAAATTTTAGTTTTGGGTCCAGAAGAACACCAAGATCATCAACTCG
etc.

OPTIONS

-i --input (default = "alignments") #Name of input folder name
-0 --output (default = "probes.txt")
-m --map.txt_file (default = "map.txt")
