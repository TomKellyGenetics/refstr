#load input data
args = commandArgs(trailingOnly=TRUE)
reference <- readLines(as.character(args[7]), warn=F)
vcf_cols <- read.table(as.character(args[6]), header=T, stringsAsFactors=F)

#split reference string into character vector
reference <- strsplit(reference, split="")[[1]]
#positions of reference
names(reference) <- 1:length(reference)

#calculate substitutions/indels - assumes all inputs are TCGG or "."
length_reference <- nchar(vcf_cols$REF)
length_alternate <- nchar(vcf_cols$ALT)
#account for missing points as "."
#length_reference[grep("[.]", vcf_cols$REF)] <- 0
length_alternate[grep("[.]", vcf_cols$ALT)] <- 0

frameshifts <- length_alternate - length_reference
adj_gtf_pos <- cumsum(frameshifts)
#table(frameshifts) # do we want to provide a text summary of frameshits detected?
# % subs + threshold of frameshifts could be warnings for new assembly / variant calling

#initialise new sequence
new_sequence <- reference
#substitute new variants into fasta sequence
for(ii in 1:nrow(vcf_cols)){
  new_sequence[vcf_cols$POS[ii]] <- vcf_cols$ALT[ii]
  if(length_reference[ii]>1) new_sequence[(vcf_cols$POS[ii]+1):(vcf_cols$POS[ii]+length_reference[ii]-1)] <- NA
  print(paste(paste(reference[vcf_cols$POS[ii]:(vcf_cols$POS[ii]+length_reference[ii]-1)], collapse=""), "substituted for", vcf_cols$ALT[ii]))
}

#test output
#new_sequence
#table(sapply(new_sequence, nchar))
#table(new_sequence)
#length(new_sequence)
#length(na.omit(new_sequence))
#table(is.na(new_sequence))

#remove empty spaces
new_sequence[grep("[.]", vcf_cols$ALT)] <- NA
new_sequence <- na.omit(new_sequence) # is this needed
#concatenate into one string
new_sequence <- paste(new_sequence, collapse='')

#test fasta
length(new_sequence)
nchar(new_sequence)

#args[4] #output name from bash
writeLines(new_sequence, con = paste0(args[4], "_string.txt"))

