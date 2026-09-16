# PA_coevolution
Python scripts to analyze raw vcf files of population level sequencing from PAO1 coevolution with OMKO1

Publication:
Kortright KE, Chan BK, Evans BR, Turner PE. Arms race and fluctuating selection dynamics in Pseudomonas aeruginosa bacteria coevolving with phage OMKO1. J Evol Biol. 2022 Nov;35(11):1475-1487. doi: 10.1111/jeb.14095. Epub 2022 Sep 27. PMID: 36168737.

Data - WGS data 2x150bp Illumina prepared with NexteraXT Library Prep Kit. Mapped to PAO1 reference (AE004091.2) with Bowtie2. Variants called with GATK.

VCFparse.py - Takes raw variant file (.vcf), filters out low-quality variants, annotates with gene information, and outputs in user friendly tab delimited [chromosome, position, reference bp, reference frequency, alternate bp, alternate frequency, gene number, gene name, quality score] for further downstream analysis.
