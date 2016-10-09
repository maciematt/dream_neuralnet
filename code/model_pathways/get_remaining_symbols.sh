#!/bin/bash


mkdir ../../data/symbols_ids

for entry in $( cat ../../data/remaining_ids.txt )
do
    ( curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id='$entry'&retmode=text' 2> /dev/null | awk -v entry=$entry 'NR == 2 {printf("%s\t%s\n", entry, $2)}' > ../../data/symbols_ids/$entry.txt )&
done


a=1
while (( a > 0 ))
do
    a=$( ls -ef | grep curl | grep -v grep | awk 'END {print NR}' )
    sleep 5
done

cat ../../data/symbols_ids/* ../../data/remaining_ids_symbols.txt


