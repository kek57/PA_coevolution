#!/usr/bin/python


import sys
import os
import glob

#os.chdir('/Users/kkortright/Desktop/Desktop/TurnerLab/2021/2021.01-January/2021.01.26-rotation_concatenated_GATK/')

genomedictfile = open('genomedict.txt','r')
genomedict = {}
for line in genomedictfile:
    line = line.strip('\n').split('\t')
    genomedict[line[0]] = [line[1], line[2], line[3]]

genomedictfile.close()

def VCFParse(argv1, argv2):
    infile = open(argv1,'r')
    outfile = open(argv2,'w')

    for line in infile:
        if line.startswith('#'):
            pass
        else:
            line = line.strip('\n').split('\t')
            chrom = line[0]
            pos = line[1]
            refbp = line[3]
            altbp = line[4]
            qual = line[5]
            low = line[6]
            test = line[8].split(':')
            if len(test) >= 2 and test[1] == "AD" and test[2] == "DP":
                dp = line[9].split(':')[2]      #read depth for alt allele
                ad = line[9].split(':')[1]      #allelic depth for ref then alt
                ref = int(ad.split(',')[0])
                alt = int(ad.split(',')[1])
                freqref = ref / (ref + alt)     #frequency of WT allele 
                freqalt = alt / (ref + alt)     #frequency of poly allele
                if pos in genomedict.keys():
                    genenumber = genomedict[pos][1]
                    genename = genomedict[pos][2]
                else:
                    genenumber = "intergenic"
                    genename = "intergenic"
                if low == "LowQual":
                    pass
                else:
                    outfile.write(str(chrom) + '\t' + str(pos) + '\t' + str(refbp) + '\t' + str(freqref) + '\t' + str(altbp) + '\t' + str(freqalt) + '\t' + str(genenumber)  + '\t' + str(genename) + '\t' + str(qual) + '\n')
            else:
                pass
    print(str('Done with ' + str(argv1[:-4]) + '!'))
    infile.close()
    outfile.close()

file1 = 'bacteria_0_raw_variants.vcf'
out1 = 'bacteria_0_parsed.txt'
VCFParse(file1, out1)

file2 = 'bacteria_PA01_raw_variants.vcf'
out2 = 'bacteria_PA01__parsed.txt'
VCFParse(file2, out2)

file3 = 'bacteria_pop1.10_raw_variants.vcf'
out3 = 'bacteria_pop1.10_parsed.txt'
VCFParse(file3, out3)

file4 = 'bacteria_pop1.1_raw_variants.vcf'
out4 = 'bacteria_pop1.1_parsed.txt'
VCFParse(file4, out4)

file5 = 'bacteria_pop1.2_raw_variants.vcf'
out5 = 'bacteria_pop1.2_parsed.txt'
VCFParse(file5, out5)

file6 = 'bacteria_pop1.3_raw_variants.vcf'
out6 = 'bacteria_pop1.3_parsed.txt'
VCFParse(file6, out6)

file7 = 'bacteria_pop1.4_raw_variants.vcf'
out7 = 'bacteria_pop1.4_parsed.txt'
VCFParse(file7, out7)

file8 = 'bacteria_pop1.5_raw_variants.vcf'
out8 = 'bacteria_pop1.5_parsed.txt'
VCFParse(file8, out8)

file9 = 'bacteria_pop1.6_raw_variants.vcf'
out9 = 'bacteria_pop1.6_parsed.txt'
VCFParse(file9, out9)

file10 = 'bacteria_pop1.7_raw_variants.vcf'
out10 = 'bacteria_pop1.7_parsed.txt'
VCFParse(file10, out10)

file11 = 'bacteria_pop1.8_raw_variants.vcf'
out11 = 'bacteria_pop1.8_parsed.txt'
VCFParse(file11, out11)

file12 = 'bacteria_pop1.9_raw_variants.vcf'
out12 = 'bacteria_pop1.9_parsed.txt'
VCFParse(file12, out12)

file13 = 'bacteria_pop2.10_raw_variants.vcf'
out13 = 'bacteria_pop2.10_parsed.txt'
VCFParse(file13, out13)

file14 = 'bacteria_pop2.1_raw_variants.vcf'
file14 = 'bacteria_pop2.1_parsed.txt'
VCFParse(file14, out14)

file15 = 'bacteria_pop2.2_raw_variants.vcf'
file15 = 'bacteria_pop2.2_parsed.txt'
VCFParse(file15, out15)

file16 = 'bacteria_pop2.3_raw_variants.vcf'
file16 = 'bacteria_pop2.3_parsed.txt'
VCFParse(file16, out16)

file17 = 'bacteria_pop2.4_raw_variants.vcf'
file17 = 'bacteria_pop2.4_parsed.txt'
VCFParse(file17, out17)

file18 = 'bacteria_pop2.5_raw_variants.vcf'
file18 = 'bacteria_pop2.5_parsed.txt'
VCFParse(file18, out18)

file19 = 'bacteria_pop2.6_raw_variants.vcf'
file19 = 'bacteria_pop2.6_parsed.txt'
VCFParse(file19, out19)

file20 = 'bacteria_pop2.7_raw_variants.vcf'
file20 = 'bacteria_pop2.7_parsed.txt'
VCFParse(file20, out20)

file21 = 'bacteria_pop2.8_raw_variants.vcf'
file21 = 'bacteria_pop2.8_parsed.txt'
VCFParse(file21, out21)

file22 = 'bacteria_pop2.9_raw_variants.vcf'
file22 = 'bacteria_pop2.9_parsed.txt'
VCFParse(file22, out22)

file23 = 'bacteria_pop3.10_raw_variants.vcf'
file23 = 'bacteria_pop3.10_parsed.txt'
VCFParse(file23, out23)

file24 = 'bacteria_pop3.1_raw_variants.vcf'
file24 = 'bacteria_pop3.1_parsed.txt'
VCFParse(file24, out24)

file25 = 'bacteria_pop3.2_raw_variants.vcf'
file25 = 'bacteria_pop3.2_parsed.txt'
VCFParse(file25, out25)

file26 = 'bacteria_pop3.3_raw_variants.vcf'
file26 = 'bacteria_pop3.3_parsed.txt'
VCFParse(file26, out26)

file27 = 'bacteria_pop3.4_raw_variants.vcf'
file27 = 'bacteria_pop3.4_parsed.txt'
VCFParse(file27, out27)

file28 = 'bacteria_pop3.5_raw_variants.vcf'
file28 = 'bacteria_pop3.5_parsed.txt'
VCFParse(file28, out28)

file29 = 'bacteria_pop3.6_raw_variants.vcf'
out29 = 'bacteria_pop3.6_parsed.txt'
VCFParse(file29, out29)

file30 = 'bacteria_pop3.7_raw_variants.vcf'
out30 = 'bacteria_pop3.7_parsed.txt'
VCFParse(file30, out30)

file31 = 'bacteria_pop3.8_raw_variants.vcf'
out31 = 'bacteria_pop3.8_parsed.txt'
VCFParse(file31, out31)

file32 = 'bacteria_pop3.9_raw_variants.vcf'
out32 = 'bacteria_pop3.9_parsed.txt'
VCFParse(file32, out32)

file33 = 'bacteria_popc.10_raw_variants.vcf'
out33 = 'bacteria_popc.10_parsed.txt'
VCFParse(file33, out33)

file34 =  'bacteria_popc.1_raw_variants.vcf'
out34 =  'bacteria_popc.1_parsed.txt'
VCFParse(file34, out34)

file35 =  'bacteria_popc.2_raw_variants.vcf'
out35 =  'bacteria_popc.2_parsed.txt'
VCFParse(file35, out35)

file36 =  'bacteria_popc.3_raw_variants.vcf'
out36 =  'bacteria_popc.3_parsed.txt'
VCFParse(file36, out36)

file37 =  'bacteria_popc.4_raw_variants.vcf'
out37 =  'bacteria_popc.4_parsed.txt'
VCFParse(file37, out37)

file38 =  'bacteria_popc.5_raw_variants.vcf'
out38 =  'bacteria_popc.5_parsed.txt'
VCFParse(file38, out38)

file39 =  'bacteria_popc.6_raw_variants.vcf'
out39 =  'bacteria_popc.6_parsed.txt'
VCFParse(file39, out39)

file40 =  'bacteria_popc.7_raw_variants.vcf'
out40 =  'bacteria_popc.7_parsed.txt'
VCFParse(file40, out40)

file41 =  'bacteria_popc.8_raw_variants.vcf'
out41 =  'bacteria_popc.8_parsed.txt'
VCFParse(file41, out41)

file42 = 'bacteria_popc.9_raw_variants.vcf'
out42 = 'bacteria_popc.9_parsed.txt'
VCFParse(file42, out42)
 


        
        
