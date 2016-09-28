*Sam playing around in bash


>get VCF file and extract info we care about

```bash
grep -E '#CHROM' -A 99999999 output1.vcf | awk 'print{$1 $2 $4 $5}' > $temp_out.list

```

>file $temp_out.list will contain the info we need for .vcf parsing, and the -A # command can be changed, i have just added a ridiculous number to ensure that it actually will get to the end of the file

> Next step is getting the .fasta file ready

```bash
sed '/>/d' ./reference.fasta | tr -d '\n' > reference_1line_noheader.fasta

```

>The above just simply removes the fasta header and then makes sure the sequence is all on one line (makes life easier later).

>appending the fasta file with changes using sed

```bash
sed '1s/\(^.\{X\}\).\{Y\}\(.*\)/\Z\2/' reference_1line_noheader.fasta

```

>Okay this will take some explaination -> in the above code X = position in fasta file (if the above steps are correct then you will only need to keep the same numbering): Y = how many positions are getting changed (in most cases it will be 1 or 2): and Z = what you are changing it to. ie:

```bash
sed '1s/\(^.\{599\}\).\{2\}\(.*\)/\CC\2/' reference_1line_noheader.fasta

```

>Would change position 599 and 600 from whatever it was to a CC
