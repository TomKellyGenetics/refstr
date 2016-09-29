#load input data
import sys
args = sys.argv

reference = open(args[7], 'r').readline()

#import numpy
#vcf_cols = numpy.genfromtxt(fname=args[6], delimiter='\t')

import numpy as np
import pandas as pd
vcf_cols = pd.read_csv(args[6], sep='\t', header=0)

#split reference string into character vector
reference = list(reference)
#positions of reference
original_position = range(1, len(reference)+1)

#calculate substitutions/indels - assumes all inputs are TCGG or "."
length_reference = []
length_alternate = []
for row in range(len(vcf_cols.index)):
     length_reference.append(len(vcf_cols.loc[row,'REF']));
     length_alternate.append(len(vcf_cols.loc[row,'ALT']));

#account for missing points as "."
for row in range(len(vcf_cols.index)):
#     if vcf_cols.loc[row,'REF'] == ".":
#             length_reference[row]=0
     if vcf_cols.loc[row,'ALT'] == ".":
             length_alternate[row]=0

#initialise new sequence
new_sequence = reference
#substitute new variants into fasta sequence
for ii in range(len(vcf_cols.index)):
    new_sequence[vcf_cols.loc[ii,'POS']-1] = vcf_cols.loc[ii,'ALT']
    if length_reference[ii] > 2:
        new_sequence[vcf_cols.loc[ii,'POS']:(vcf_cols.loc[ii,'POS']+length_reference[ii]-1)] = "NA"
    elif length_reference[ii] > 1:
        new_sequence[vcf_cols.loc[ii,'POS']] = "NA"
    print vcf_cols.loc[ii,'REF'],"substituted for",vcf_cols.loc[ii,'ALT']

#test output
#new_sequence
#table(sapply(new_sequence, len))
#table(new_sequence)
#length(new_sequence)
#length(na.omit(new_sequence))
#table(is.na(new_sequence))

#remove empty spaces
new = []
for seq in new_sequence:
    if seq != "NA":
        new.append(seq)
new_sequence = new
#concatenate into one string
new_sequence = ''.join(new_sequence)

#test fasta
len(new_sequence)

#args[4] #output name from bash
open(''.join([args[4], "_string.txt"]), 'w+').write(new_sequence)

#read gtf data
#gtf <- data.table::fread(as.character(args[3]), data.table = F)

#calculate framshifts
#frameshifts <- length_alternate - length_reference
#adj_gtf_pos <- cumsum(frameshifts[order(vcf_cols.loc[:,'POS'])])
#table(frameshifts) # do we want to provide a text summary of frameshits detected?
# % subs + threshold of frameshifts could be warnings for new assembly / variant calling

#original fasta pos
#original_position
#vcf_cols.loc[:,'POS']

#add starting point to new positions
#if(vcf_cols$POS[1]!=1){
#    position <- c(1, vcf_cols.loc[:,'POS'][order(vcf_cols.loc[:,'POS'])])
#    adj_gtf_pos <- c(0, adj_gtf_pos)
#  } else {
#   position <- vcf_cols.loc[:,'POS'][order(vcf_cols$POS)]
#}
#if(position[length(position)]!=length(reference)){
#  position <- c(position, length(reference))
#} 


#new pos
#new_positions <- as.numeric(original_position)
#for(ii in 1:(length(position)-1)){
#  new_positions[position[ii]:(position[ii+1]-1)] <-new_positions[position[ii]:(position[ii+1]-1)] + adj_gtf_pos[ii]
#}
#new_positions[length(new_positions)] <- len(new_sequence)

#substitute new gtf positions #do we want to warn if frameshift within a gene?
#new_gtf <- gtf
#new_gtf$V4 <- new_positions[gtf$V4]
#new_gtf$V5 <- new_positions[gtf$V5]

#output gtf file
#write.table(new_gtf, file = paste0(args[4], ".gtf"), header=F, quotes=F)

