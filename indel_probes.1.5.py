#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This program finds oligonucleotide probe sequencing suitable for 
detecting polymorphisms using gDNA hybridization to probes on glass slides.
For each position, one probe is made to straddle the deletion and one is 
specific to the insertion.
s1ple screens are performed to exclude runs of NNs and s1ple repeat regions.

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
-min --min_indel (default = 10) #Minimum length of insertion/deletion events
-max --max_indel (default = 100) #Maximum length of insertion/deletion events
-GC --minGCcontect (default = 35%)


v.1.0 
I imposed a fairly stringent screen against bad alignments by excluding any
sequences with a second deletion within 50bp on either side and by
excluding any sequences with N's
sequences with di and trinucleotide repeats are removed
v.1.1 probe names now reflect chromosome and genome position
v. 1.5 updated to Bio SeqIO, 
"""


from string import *
from Bio import SeqIO
from oligo_functions_1_1 import *
import argparse
import fileinput

parser = argparse.ArgumentParser(description='Make oligo probes from aligned gDNA sequences')
parser.add_argument('-i', default = 'alignments', help='input folder name',type=str)
parser.add_argument('-o', default = 'probes.txt', help='output file name',type=argparse.FileType('w'))
parser.add_argument('-m', default = 'map.txt', help='output file name',type=argparse.FileType('r'))
#parser.add_argument('-min', default = 10, help='minimum length of in/del event',type=int)
#parser.add_argument('-max', default = 100, help='maximum length of in/del event',type=int)
#parser.add_argument('-GC', default = 35, help='minimum proportion of GC in 100bp region surrounding indel',type=int)

args = parser.parse_args()
input_folder = args.i
output_file = args.o
mapping = args.m


#parser.print_help() 

#output_file = open('probes.txt','w')
#output_file_fasta = open('probes.fasta.txt','w')
#output_probe_translation = open('probe_position_translation.txt','w') #I used the wrong genomic map position in original names for array, forgot to delete deletions. Use this translation table to replace names with correct names
#mapping = open('map.txt','r')    

in_seq_oligo_name = []
del_seq_oligo_name = []
in_seq_oligo = []
del_seq_oligo = []


#folder = 1
for line in mapping:
#while folder < 20:
    #get genome position of folder sequence
#    line = mapping.readline()
    data = line.split('\t')
    file_name = data[0]
    chromosome = data[1]
    chromosome_number = chromosome[3:5]
    genome_position = int(data[2])
    
    handle = open(input_folder + '/' + str(file_name), 'rU')
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    #read s1 and s2h sequences
#    x = 0
#    for cur_record in SeqIO.parse(input_file, "fasta"):
        
#        if x == 0:
 #           s1_region_name = cur_record.title.split(' ')[0]
    s1name = records[0].description
    s1_seq_name= s1name.split()[0]
    s2name = records[1].description
    s2_seq_name= s2name.split()[0]
    s1 = records[0].seq
    s1_seq = str(s1)
  #          x += 1
   #     else:
    #        s2_region_name = cur_record.title.split(' ')[0]
    s2 = records[1].seq
    s2_seq = str(s2)

    #first do s1 deletions
    start_position = 0
    stop_position = 0
    count = 0
    n = 0
    while n < len(s1_seq):
        if s1_seq[n] == '-': #determine the length of the deletion in s1 sequence
            count += 1
            if count == 1:
                start_position = n
        else:
            if 9 < count < 100: #if the in/del is between 10 and 100 bp long
                stop_position = n                
                s1_region = s1_seq[start_position-50:stop_position+50]
                number_s1_gaps = s1_region.count('-')
                s2_region = s2_seq[start_position-50:stop_position+50]
                number_s2_gaps = s2_region.count('-')
                number_GC = s1_region.count('G') + s1_region.count('C')
                if number_s1_gaps > count:pass
                elif number_GC<35:pass
                elif number_s2_gaps > 0:pass
                elif start_position < 50:pass
                elif len(s1_seq)-stop_position < 50:pass
                elif 'N' in s1_region:pass
                elif 'ATATATATAT' in s1_region:pass
                elif 'ACACACACAC' in s1_region:pass
                elif 'AGAGAGAGAG' in s1_region:pass
                elif 'CGCGCGCGCG' in s1_region:pass
                elif 'CTCTCTCTCT' in s1_region:pass
                elif 'GTGTGTGTGT' in s1_region:pass
                elif 'AATAATAAT' in s1_region:pass
                elif 'ATTATTATT' in s1_region:pass
                elif 'AACAACAAC' in s1_region:pass
                elif 'ACCACCACC' in s1_region:pass
                elif 'AAGAAGAAG' in s1_region:pass
                elif 'AGGAGGAGG' in s1_region:pass
                elif 'TTCTTCTTC' in s1_region:pass
                elif 'TCCTCCTCC' in s1_region:pass
                elif 'TTGTTGTTG' in s1_region:pass
                elif 'TGGTGGTGG' in s1_region:pass
                elif 'CCGCCGCCG' in s1_region:pass
                elif 'CGGCGGCGG' in s1_region:pass
                elif 'ATCATCATC' in s1_region:pass
                elif 'ATGATGATG' in s1_region:pass
                elif 'TACTACTAC' in s1_region:pass
                elif 'TAGTAGTAG' in s1_region:pass
                elif 'TCGTCGTCG' in s1_region:pass
                elif 'TGCTGCTGC' in s1_region:pass
                elif 'GTCGTCGTC' in s1_region:pass
                elif 'CGACGACGA' in s1_region:pass
                elif 'CAGCAGCAG' in s1_region:pass
                elif 'AAATAAAT' in s1_region:pass
                elif 'AATTAATT' in s1_region:pass
                elif 'AAAAAAAA' in s1_region:pass
                elif 'TTTTTTTT' in s1_region:pass
                elif 'N' in s2_region:pass
                elif 'ATATATATAT' in s2_region:pass
                elif 'ACACACACAC' in s2_region:pass
                elif 'AGAGAGAGAG' in s2_region:pass
                elif 'CGCGCGCGCG' in s2_region:pass
                elif 'CTCTCTCTCT' in s2_region:pass
                elif 'GTGTGTGTGT' in s2_region:pass
                elif 'AATAATAAT' in s2_region:pass
                elif 'ATTATTATT' in s2_region:pass
                elif 'AACAACAAC' in s2_region:pass
                elif 'ACCACCACC' in s2_region:pass
                elif 'AAGAAGAAG' in s2_region:pass
                elif 'AGGAGGAGG' in s2_region:pass
                elif 'TTCTTCTTC' in s2_region:pass
                elif 'TCCTCCTCC' in s2_region:pass
                elif 'TTGTTGTTG' in s2_region:pass
                elif 'TGGTGGTGG' in s2_region:pass
                elif 'CCGCCGCCG' in s2_region:pass
                elif 'CGGCGGCGG' in s2_region:pass
                elif 'ATCATCATC' in s2_region:pass
                elif 'ATGATGATG' in s2_region:pass
                elif 'TACTACTAC' in s2_region:pass
                elif 'TAGTAGTAG' in s2_region:pass
                elif 'TCGTCGTCG' in s2_region:pass
                elif 'TGCTGCTGC' in s2_region:pass
                elif 'GTCGTCGTC' in s2_region:pass
                elif 'CGACGACGA' in s2_region:pass
                elif 'CAGCAGCAG' in s2_region:pass
                elif 'AAATAAAT' in s2_region:pass
                elif 'AATTAATT' in s2_region:pass
                elif 'AAAAAAAA' in s2_region:pass
                elif 'TTTTTTTT' in s2_region:pass
                else:
#                    print folder
#                    print "s1_del"
#                    print s1_region
#                    print s2_region
                    in_seq_oligo = get_in_oligo(s2_region)
 #                   print "mark"
 #                   print in_seq_oligo
                    del_seq_oligo = get_del_oligo(s1_region)
                                        
                    s1_seq_from_beginning = s1_seq[0:start_position] 
                    total_number_s1_gaps = s1_seq_from_beginning.count('-')
                    position_in_contig = start_position - 1 - total_number_s1_gaps
                    
#                    in_seq_oligo_name = get_in_oligo_name(chromosome_number,genome_position+start_position,'s2')
#                    del_seq_oligo_name = get_del_oligo_name(chromosome_number,genome_position+start_position,'s1')
    
                    real_del_oligo_position = get_del_oligo_name(chromosome_number,genome_position+position_in_contig,s1_seq_name)
                    real_in_oligo_position = get_in_oligo_name(chromosome_number,genome_position+position_in_contig,s2_seq_name)
                    
#                    print in_seq_oligo_name
#                    print in_seq_oligo
#                    print del_seq_oligo_name
#                    print del_seq_oligo
                    
                    if del_seq_oligo == "skip": pass
                    elif in_seq_oligo == "skip": pass
                    else:
                        output_file.write('%s\t%s\n' %(real_del_oligo_position,del_seq_oligo))
                        output_file.write('%s\t%s\n' %(real_in_oligo_position,in_seq_oligo))
                        
#                        output_file_fasta.write('>%s\n%s\n' %(real_del_oligo_position,s1_seq[start_position-300:stop_position+300]))
#                        output_file_fasta.write('>%s\n%s\n' %(real_in_oligo_position,s2_seq[start_position-300:stop_position+300]))
                    
#                        output_probe_translation.write('%s\t%s\t%s\n' %(del_seq_oligo_name,real_del_oligo_position,del_seq_oligo))
#                        output_probe_translation.write('%s\t%s\t%s\n' %(in_seq_oligo_name,real_in_oligo_position,in_seq_oligo))


            count = 0
            start_position = 0
            stop_position = 0
        n += 1
#    print n
#    print len(s1_seq)
    #then do s2h sequences
    start_position = 0
    stop_position = 0
    count = 0
    n = 0
    while n < len(s2_seq):
        if s2_seq[n] == '-':
            count += 1
            if count == 1:
                start_position = n
        else:
            if 9 < count < 100: #if the in/del is between 10 and 100 bp long 
                stop_position = n
                s1_region = s1_seq[start_position-50:stop_position+50]
                s2_region = s2_seq[start_position-50:stop_position+50]
                number_s1_gaps = s1_region.count('-')
                number_s2_gaps = s2_region.count('-')
                number_GC = s2_region.count('G') + s2_region.count('C')
                if number_s2_gaps > count:pass
                elif number_GC<35:pass
                elif number_s1_gaps > 0:pass
                elif start_position < 50:pass
                elif len(s1_seq)-stop_position < 50:pass
                elif 'N' in s1_region:pass
                elif 'ATATATATAT' in s1_region:pass
                elif 'ACACACACAC' in s1_region:pass
                elif 'AGAGAGAGAG' in s1_region:pass
                elif 'CGCGCGCGCG' in s1_region:pass
                elif 'CTCTCTCTCT' in s1_region:pass
                elif 'GTGTGTGTGT' in s1_region:pass
                elif 'AATAATAAT' in s1_region:pass
                elif 'ATTATTATT' in s1_region:pass
                elif 'AACAACAAC' in s1_region:pass
                elif 'ACCACCACC' in s1_region:pass
                elif 'AAGAAGAAG' in s1_region:pass
                elif 'AGGAGGAGG' in s1_region:pass
                elif 'TTCTTCTTC' in s1_region:pass
                elif 'TCCTCCTCC' in s1_region:pass
                elif 'TTGTTGTTG' in s1_region:pass
                elif 'TGGTGGTGG' in s1_region:pass
                elif 'CCGCCGCCG' in s1_region:pass
                elif 'CGGCGGCGG' in s1_region:pass
                elif 'ATCATCATC' in s1_region:pass
                elif 'ATGATGATG' in s1_region:pass
                elif 'TACTACTAC' in s1_region:pass
                elif 'TAGTAGTAG' in s1_region:pass
                elif 'TCGTCGTCG' in s1_region:pass
                elif 'TGCTGCTGC' in s1_region:pass
                elif 'GTCGTCGTC' in s1_region:pass
                elif 'CGACGACGA' in s1_region:pass
                elif 'CAGCAGCAG' in s1_region:pass
                elif 'AAATAAAT' in s1_region:pass
                elif 'AATTAATT' in s1_region:pass
                elif 'AAAAAAAA' in s1_region:pass
                elif 'TTTTTTTT' in s1_region:pass
                elif 'N' in s2_region:pass
                elif 'ATATATATAT' in s2_region:pass
                elif 'ACACACACAC' in s2_region:pass
                elif 'AGAGAGAGAG' in s2_region:pass
                elif 'CGCGCGCGCG' in s2_region:pass
                elif 'CTCTCTCTCT' in s2_region:pass
                elif 'GTGTGTGTGT' in s2_region:pass
                elif 'AATAATAAT' in s2_region:pass
                elif 'ATTATTATT' in s2_region:pass
                elif 'AACAACAAC' in s2_region:pass
                elif 'ACCACCACC' in s2_region:pass
                elif 'AAGAAGAAG' in s2_region:pass
                elif 'AGGAGGAGG' in s2_region:pass
                elif 'TTCTTCTTC' in s2_region:pass
                elif 'TCCTCCTCC' in s2_region:pass
                elif 'TTGTTGTTG' in s2_region:pass
                elif 'TGGTGGTGG' in s2_region:pass
                elif 'CCGCCGCCG' in s2_region:pass
                elif 'CGGCGGCGG' in s2_region:pass
                elif 'ATCATCATC' in s2_region:pass
                elif 'ATGATGATG' in s2_region:pass
                elif 'TACTACTAC' in s2_region:pass
                elif 'TAGTAGTAG' in s2_region:pass
                elif 'TCGTCGTCG' in s2_region:pass
                elif 'TGCTGCTGC' in s2_region:pass
                elif 'GTCGTCGTC' in s2_region:pass
                elif 'CGACGACGA' in s2_region:pass
                elif 'CAGCAGCAG' in s2_region:pass
                elif 'AAATAAAT' in s2_region:pass
                elif 'AATTAATT' in s2_region:pass
                elif 'AAAAAAAA' in s2_region:pass
                elif 'TTTTTTTT' in s2_region:pass
                else:
#                    print folder
#                    print "s2_del"
#                    print s1_region
#                    print s2_region

                    in_seq_oligo = get_in_oligo(s1_region)
                    del_seq_oligo = get_del_oligo(s2_region)

                    s1_seq_from_beginning = s1_seq[0:start_position]
                    total_number_s1_gaps = s1_seq_from_beginning.count('-')
                    position_in_contig = start_position - total_number_s1_gaps

#                    in_seq_oligo_name = get_in_oligo_name(chromosome_number,genome_position+start_position,'s1')
 #                   del_seq_oligo_name = get_del_oligo_name(chromosome_number,genome_position+start_position,'s2')

                    real_in_oligo_position = get_in_oligo_name(chromosome_number,genome_position+position_in_contig,s1_seq_name)
                    real_del_oligo_position = get_del_oligo_name(chromosome_number,genome_position+position_in_contig,s2_seq_name)

                    
#                    print in_seq_oligo_name
#                    print in_seq_oligo
#                    print del_seq_oligo_name
#                    print del_seq_oligo

                    if del_seq_oligo == "skip": pass
                    elif in_seq_oligo == "skip": pass
                    else:

                        output_file.write('%s\t%s\n' %(real_in_oligo_position,in_seq_oligo))
                        output_file.write('%s\t%s\n' %(real_del_oligo_position,del_seq_oligo))
                        
#                        output_file_fasta.write('>%s\n%s\n' %(real_in_oligo_position,s1_seq[start_position-300:stop_position+300]))
#                        output_file_fasta.write('>%s\n%s\n' %(real_del_oligo_position,s2_seq[start_position-300:stop_position+300]))

#                        output_probe_translation.write('%s\t%s\t%s\n' %(in_seq_oligo_name,real_in_oligo_position,in_seq_oligo))
#                        output_probe_translation.write('%s\t%s\t%s\n' %(del_seq_oligo_name,real_del_oligo_position,del_seq_oligo))

            count = 0
            start_position = 0
            stop_position = 0
        n += 1
#    print n
#    print len(s1_seq)
    print file_name
#    folder += 1


#    input_file.close()
    
output_file.close()
#output_file_fasta.close()
#output_probe_translation.close()
mapping.close()
