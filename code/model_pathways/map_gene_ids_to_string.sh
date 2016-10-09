#!/bin/bash


list_of_genes=$( cat $1 )

human_species_id=9606

curl -d "identifiers="$list_of_genes"\&amp;species="$human_species_id string-db.org/api/tsv/resolveList > $2

