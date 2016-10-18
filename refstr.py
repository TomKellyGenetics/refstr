#load input data
import sys
args = sys.argv

reference = open(args[7], 'r').readline()

#import numpy
#vcf_cols = numpy.genfromtxt(fname=args[6], delimiter='\t')

import numpy as np
import pandas as pd
vcf_cols = pd.read_csv(args[6], sep='\t', header=0)
vcf_cols = vcf_cols.sort()

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
    new_sequence[(vcf_cols.loc[ii,'POS']-1)] = vcf_cols.loc[ii,'ALT']
    if length_reference[ii] > 2:
        new_sequence[vcf_cols.loc[ii,'POS']:(vcf_cols.loc[ii,'POS']+length_reference[ii]-1)] = None
    elif length_reference[ii] > 1:
        new_sequence[vcf_cols.loc[ii,'POS']] = None
    print vcf_cols.loc[ii,'REF'],"substituted for",vcf_cols.loc[ii,'ALT']

#test output
#new_sequence
#table(sapply(new_sequence, len))
#table(new_sequence)
#length(new_sequence)
#length(na.omit(new_sequence))
#table(is.na(new_sequence))


#concatenate into one string
new_sequence = ''.join(new_sequence)
#remove empty spaces
new_sequence.replace("n", "")
#new = []
#for seq in new_sequence:
#    if seq != "None":
#        new.append(seq)
#new_sequence = new

#test fasta
len(new_sequence)

#args[4] #output name from bash
open(''.join([args[4], "_string.txt"]), 'w+').write(new_sequence)

#read gtf data
gtf = pd.read_csv(args[3], sep='\t', header=0)

#calculate frameshifts
frameshifts = np.array(length_alternate) - np.array(length_reference)
adj_gtf_pos = frameshifts.cumsum()
#table(frameshifts) # do we want to provide a text summary of frameshits detected?
# % subs + threshold of frameshifts could be warnings for new assembly / variant calling

#original fasta pos
#original_position
#vcf_cols.loc[:,'POS']

#add starting point to new positions
position = vcf_cols.loc[:,'POS']
if vcf_cols[:POS][1] != 1:
    position = pd.Series([1])
    position.append(vcf_cols.loc[:,'POS'])
    temp = [1]
    temp.extend(adj_gtf_pos)
    adj_gtf_pos = temp
    temp = [1]
    temp.extend(frameshifts)
    frameshifts = temp
position = position.append(pd.Series(len(reference)))


#new pos
#new_positions = []
#for jj in range(len(position)-1):
#    if frameshifts[jj] < 0:
#        for i in range(1,-frameshifts[jj]+1)
#            new_positions.append(position[jj]+adj_gtf_pos[jj-1]) #repeats for del
#        new_positions.append(range(-frameshifts[jj]+position[jj]+adj_gtf_pos[jj]),(position[jj+1]adj_gtf_pos[jj])))
#    if frameshifts[jj] > 0:
#        new_positions.append(position[jj]+adj_gtf_pos[jj]) #skip ins
#        new_positions.append(range(frameshifts[jj]+position[jj]+adj_gtf_pos[jj-1]+1),(position[jj+1]+adj_gtf_pos[jj]))) #skip ins
#    if frameshifts[jj] == 0:
#        new_positions.append(range(position[jj]+adj_gtf_pos[jj]),(position[jj+1]+adj_gtf_pos[jj]))) #keep pos for subs
#    print(jj)
#


#new_positions = original_position
#for ii in range(len(position)-1):


new_positions = []
for ii in range(len(position)-1):
    new_positions.append(map(lambda x: x+adj_gtf_pos[ii],  original_position[position[ii]:(position[ii+1]-1)]))
new_positions.append(range(position[len(position)-1], len(new_sequence)+1)+adj_gtf_pos[len(position)-1])
new_positions = [item for sublist in new_positions for item in sublist]

#substitute new gtf positions #do we want to warn if frameshift within a gene?
new_gtf = gtf
new_gtf[:, 3] = new_positions[gtf[:, 3]]
new_gtf[:, 4] = new_positions[gtf[:, 4]]

#output gtf file
write.table(new_gtf, file = paste0(args[4], ".gtf"), header=F, quotes=F)

