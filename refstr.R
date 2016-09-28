#load inpute data
reference <- data.table::fread("reference_string.txt", data.table = F)
vcf_cols <- data.table::fread("output.txt", data.table = F)

#split reference string into character vector
reference <- strsplit(reference, split="")[[1]]
#positions of reference
names(reference) <- 1:length(reference)



head(vcf_cols)
