#! \bin\bash
# args: input_vcf($1) ref_seq($2) ref_gtf($3) output_name($4) new_fasta_header($5) cores($6))

#echo $1 input
#echo $2 reference
#echo ${6:-1} cores (defaults to 1)

#define intermediate file names
reference_string=$(echo $2 | sed -e "s/.fasta/_string.txt/")
#echo $reference_string
input_tsv=$(echo $1 | sed -e "s/.vcf/.tsv/")
#echo $input_tsv

#tail removes header, sed removes lines breaks
#tail -n +2 $reference_string > $reference_string
#sed ':a;N;$!ba;s/\n//g' "${reference_string"
#replace with sed to remove header and tr to remove line breaks (preserves variable output)
sed '/>/d' ./${2} | tr -d '\n' > $reference_string


#remove header from VCF
grep "#CHROM" $1 > $input_tsv
sed -n -E -e '/#CHROM/,$ p' $1 | sed '1 d' >> $input_tsv
sed -i 's/#//g' $input_tsv

#Rscript refstr.R "${@:1:5}" $input_tsv $reference_string ${6:-1} #initialise R script and pass pass arguments
python refstr.py "${@:1:5}" $input_tsv $reference_string ${6:-1} #initialise R script and pass pass arguments

#cut fasta to rows of 60 chars
fold -w 60 -s ${4}_string.txt > ${4}.fasta
#add header
sed -i "1s/^/\>${5}\n/" ${4}.fasta
