#! \bin\bash
# args: input_vcf ref_seq output_name

#tail removes header, sed removes lines breaks
tail -n +2 $2 | sed ':a;N;$!ba;s/\n//g' > reference_string.txt

#remove header from VCF
grep "#CHROM" output1.vcf > output.txt
sed -n -E -e '/#CHROM/,$ p' output1.vcf| sed '1 d' >> output.txt