#load inpute data
reference <- readLines("reference_string.txt")
vcf_cols <- datatable::fread("output.txt")

#split reference string into character vector
reference <- strsplit(reference, split="")[[1]]
#positions of reference
names(reference) <- 1:length(reference)



head(vcf_cols)
