#load inpute data
reference <- readLines("reference_string.txt")
vcf_cols <- data.table::fread("output.txt", data.table = F)

#split reference string into character vector
reference <- strsplit(reference, split="")[[1]]
#positions of reference
names(reference) <- 1:length(reference)

#calculate substitutions/indels - assumes all inputs are TCGG or "."
length_reference <- nchar(vcf_cols$REF)
length_alternate <- nchar(vcf_cols$ALT)
#account for missing points as "."
length_reference [grep("[.]", vcf_cols$REF)] <- 0
length_alternate [grep("[.]", vcf_cols$ALT)] <- 0

frameshifts <- length_alternate - length_reference
table(frameshifts) # do we want to provide a text summary of frameshits detected?
# % subs + threshold of frameshifts could be warnings for new assembly / variant calling

#initialise new sequence
new_sequence <- reference
#substitute new variants into fasta sequence
for(ii in 1:nrow(vcf_cols)){
  new_sequence[vcf_cols$POS[ii]:(vcf_cols$POS[ii]+length_reference[ii]-1)] <- vcf_cols$ALT[ii]
  print(paste(reference[vcf_cols$POS[ii]:(vcf_cols$POS[ii]+length_reference[ii]-1)], "substituted for", vcf_cols$ALT[ii]))
}

head(vcf_cols)

Sys.getenv('3') #output name from bash
