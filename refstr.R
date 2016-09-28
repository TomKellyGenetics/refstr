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

frameshifts <- ength_alternate - length_reference

for(ii in 1:nrow(vcf_cols)){
  
}

head(vcf_cols)
