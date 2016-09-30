#!/usr/bin/julia

#load input data
reference = readline(open(ARGS[7]))
#Pkg.add("DataFrames")
using DataFrames
vcf_cols = readtable(open(ARGS[6]), separator='\t')
vcf_cols = sort(vcf_cols)

#split reference string into character vector
new_ref = Any[]
for ii in 1:length(reference)
    push!(new_ref, reference[ii])
    end
reference = new_ref
#positions of reference
original_position = 1:length(reference)

#calculate substitutions/indels - assumes all inputs are TCGG or "."
length_reference = Any[]
length_alternate = Any[]
for ii in 1:nrow(vcf_cols)
    push!(length_reference, length(vcf_cols[:REF][ii]))
    push!(length_alternate, length(vcf_cols[:ALT][ii]))
    end


#account for missing points as "."
for ii in 1:nrow(vcf_cols)
#     if vcf_cols[:REF][ii] == "."
#        length_reference[ii] = 0
#        end
    if vcf_cols[:ALT][ii] == "."
        length_alternate[ii] = 0
        print(ii)
        end
    end

#initialise new sequence
new_sequence = reference
#substitute new variants into fasta sequence
for ii in 1:nrow(vcf_cols)
    new_sequence[vcf_cols[:POS][ii]] = vcf_cols[:ALT][ii]
        if length_reference[ii] > 2
            new_sequence[(vcf_cols[:POS][ii]+1):(vcf_cols[:POS][ii]+length_reference[ii])] = "" #NA
        elseif length_reference[ii] > 1
            new_sequence[(vcf_cols[:POS][ii]+1)] = "" #NA
            end
    println(vcf_cols[:REF][ii]," substituted for ",vcf_cols[:ALT][ii])    
    end

#test output
#new_sequence
#table(sapply(new_sequence, nchar))
#table(new_sequence)
#length(new_sequence)
#length(na.omit(new_sequence))
#table(is.na(new_sequence))

#remove empty space
#new = Any[]
#for seq in new_sequence
#    if seq != NA
#        push!(new, seq)
#        end
#    end
#new_sequence = new
#concatenate into one string
new_sequence = join(new_sequence)

#test fasta
#length(new_sequence)
#nchar(new_sequence)

#ARGS[4] #output name from bash
touch(join([ARGS[4], "_string.txt"]))
write(open(join([ARGS[4], "_string.txt"]), "w"), new_sequence)


#read gtf data
gtf = readtable(open(ARGS[3]), separator='\t')


#calculate frameshifts
frameshifts = length_alternate - length_reference
adj_gtf_pos = cumsum(frameshifts)
#table(frameshifts) # do we want to provide a text summary of frameshits detected?
# % subs + threshold of frameshifts could be warnings for new assembly / variant calling

#original fasta pos
#original_position
#vcf_cols[:POS]


#add starting point to new positions
position = vcf_cols[:POS]
if vcf_cols[:POS][1] != 1
    unshift!(position, 1)
    unshift!(adj_gtf_pos, 0)
    unshift!(frameshifts, 0)
    end
push!(position, length(reference))

#new pos
new_positions = Any[]
for jj in 1:(length(position)-1)
           if frameshifts[jj] < 0
               for i in 1:(-frameshifts[jj])
                   push!(new_positions, position[jj]+adj_gtf_pos[jj-1]) #repeats for del
                   end
               append!(new_positions, (-frameshifts[jj]+position[jj]+adj_gtf_pos[jj]):(position[jj+1]-1+adj_gtf_pos[jj]))
               end
           if frameshifts[jj] > 0
               push!(new_positions, position[jj]+adj_gtf_pos[jj]) #skip ins
               append!(new_positions, (frameshifts[jj]+position[jj]+adj_gtf_pos[jj-1]+1):(position[jj+1]+adj_gtf_pos[jj]-1)) #skip ins
               end
           if frameshifts[jj] == 0
               append!(new_positions, (position[jj]+adj_gtf_pos[jj]):(position[jj+1]+adj_gtf_pos[jj]-1)) #keep pos for subs
               end
           print(jj)
           end

#substitute new gtf positions #do we want to warn if frameshift within a gene?
new_gtf = gtf
new_gtf[:, 4] = new_positions[gtf[:4]]
new_gtf[:, 5] = new_positions[gtf[:5]]

#output gtf file
touch(join([ARGS[4], ".gtf"]))
writetable(join([ARGS[3], ".gtf"]), new_gtf, separator='\t')
