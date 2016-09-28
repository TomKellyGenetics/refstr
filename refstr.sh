#! \bin\bash
# args: input_vcf($1) ref_seq($2) output_name($3) new_fasta_header($4)

#tail removes header, sed removes lines breaks
tail -n +2 $2 | sed ':a;N;$!ba;s/\n//g' > reference_string.txt

#remove header from VCF
grep "#CHROM" $1 > output.tsv
sed -n -E -e '/#CHROM/,$ p' $1 | sed '1 d' >> output.tsv
sed -i 's/#//g' output.tsv

Rscript refstr.R "$@" #initialise R script and pass pass arguments

#cut fasta to rows of 60 chars
fold -w 60 -s ${3}_string.txt > ${3}.fasta
#add header
echo sed -i "1s/^/\>${4}\n/" output.fasta