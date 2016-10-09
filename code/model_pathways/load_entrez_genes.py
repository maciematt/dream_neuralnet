#!/usr/bin/env python2.7



import os.path, mygene
import copy
import cPickle as pickle
from subprocess import call
#from Bio import Entrez




def load_entrez_id_to_symbol():


    entrez_id_to_symbol_file = '../../data/entrez_id_to_symbol.pkl'


    if os.path.isfile(entrez_id_to_symbol_file):
        entrez_id_to_symbol = pickle.load(open(entrez_id_to_symbol_file, 'rb'))

    else:

        location_biocarta_kegg_reactome = '../../src/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt'
        location_entrez = '../../src/PASCAL/resources/genesets/msigdb/msigdb.v4.0.entrez.gmt'


        entrez_gene_ids = set()

        for location in [location_entrez, location_biocarta_kegg_reactome]:
            with open(location, 'r') as location_f:
                for line in location_f.readlines():
                    line = line.strip()
                    l = line.split('\t')
                    entrez_gene_ids = entrez_gene_ids | set([int(x) for x in l[2:]])


        entrez_gene_ids = list(entrez_gene_ids)


        entrez_id_to_symbol = {}

        mg = mygene.MyGeneInfo()
        entrez_id_to_symbol_mg = mg.querymany(entrez_gene_ids, scopes = 'entrezgene', fields = ['entrezgene', 'symbol'], species = 'human')



        def filter_entrez(x):
            try:
                return {'gene_id' : x['entrezgene'], 'gene_symbol' : x['symbol']}
            except KeyError:
                return None



        entrez_id_to_symbol_mg = [y for y in (filter_entrez(x) for x in entrez_id_to_symbol_mg) if y is not None]


        remaining_ids = set(entrez_gene_ids) - set([x['gene_id'] for x in entrez_id_to_symbol_mg])

        with open('../../data/remaining_ids.txt', 'w') as remaining_ids_fhandle:
            remaining_ids_fhandle.write('\n'.join(remaining_ids) + '\n')

        #Entrez.email = 'maciejewski.matt@gmail.com'

        #remaining_symbols_handles = [Entrez.efetch(db = 'gene', id = x, retmode = 'text') for x in remaining_ids]
        #pickle.dump(remaining_symbols_handles, open(entrez_id_to_symbol_file, 'wb'))


        call(['./get_remaining_symbols.sh'])


        with open('../../data/remaining_ids_symbols.txt', 'r') as remaining_symbols_ids:
            remaining_symbols = {int(y[0]) : y[1] for y in [x.strip().split('\t') for x in remaining_symbols_ids.readlines()]}


        entrez_id_to_symbol = {x['gene_id'] : x['gene_symbol'] for x in entrez_id_to_symbol_mg}
        entrez_id_to_symbol.update(remaining_symbols)


        pickle.dump(entrez_id_to_symbol, open(entrez_id_to_symbol_file, 'wb'))


    return(entrez_id_to_symbol)

