#################################################################
# By: Jessica Qiu
# Date: December 15, 2017
# Purpose: To format results for trinotate
# How to Use: ./format_create_map.sh <protein>.fasta <pep>.fasta
#################################################################

#!/bin/bash

## Extracts Only Header Name
awk '{print $1}' $1 > edited.$1
awk '{print $1}' $2 > edited.$2

## Count Lengths and Summate
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' edited.$2 > pep_lengths.txt

## Loop to go through the pep_lengths.txt to get total length
count=0
counter=0
check=1
cat pep_lengths.txt | while read i
do
   line=$(($count % 2))

   if (( $count >= 2 )) && (( $line == 0 ))
   then
      name=$(echo $i)
      export ORIGINAL=$(echo "${name//>}")
      original[check]=$ORIGINAL

      NEW_HEADER=$(echo -e "${name} len:${length} (+)  ${ORIGINAL}:1-${length}(+)")
      new_header[$counter]=$NEW_HEADER

      sed "s/${original[$counter]}/${new_header[$counter]}/" edited.$2 >> new_headers.fasta
      ((counter++))
      ((check++))
   fi

   if (( $count == 0 ))
   then
      name=$(echo $i)

      ORIGINAL=$(echo "${name//>}")
      original[0]=$ORIGINAL
   fi

   if (( $count % 2 == 1 ))
   then
      export length=$(echo $i)
   fi

   ((count++));
done

## Heredoc to create the gene_trans_map
(
cat << MAP
#!/usr/bin/perl

open(FASTA, "<edited.$1");
while(<FASTA>) {
    chomp(\$_);
    if (\$_ =~  m/^>/ ) {
        my \$header = \$_;
        substr(\$header,0,1)="";
        print "\$header\t\$header\n";
    }
}
MAP
) > create_gene_trans_map.pl

perl create_gene_trans_map.pl > $1.gene_trans_map
mv "$1.gene_trans_map" "$(echo $1.gene_trans_map | sed -r 's/.fasta//')"
