#!/usr/bin/env python2.7



import os.path
import pandas as pd
import cPickle as pickle
from pprint import pprint
from load_entrez_genes import load_entrez_id_to_symbol





class Pathway:

    def __init__(self, name, link, entrez_gene_ids, entrez_id_to_symbol):

        self.name = name
        self.link = link
        self.entrez_gene_ids = [int(x) for x in entrez_gene_ids]
        self.gene_symbols = [entrez_id_to_symbol[x] for x in self.entrez_gene_ids]




class MsigDB:

    def __init__(self):

        location_msigdb = '../../data/msigdb.pkl'

        self.entrez_id_to_symbol = load_entrez_id_to_symbol()

        if os.path.isfile(location_msigdb):
            temp_dict = pickle.load(open(location_msigdb, 'rb'))
            self.__dict__.update(temp_dict)

        else:
            self.location_biocarta_kegg_reactome = '../../src/PASCAL/resources/genesets/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt'
            self.location_entrez = '../../src/PASCAL/resources/genesets/msigdb/msigdb.v4.0.entrez.gmt'

            self.pathways_entrez = self.slurp_data(self.location_entrez)
            self.pathways_biocarta_kegg_reactome = self.slurp_data(self.location_biocarta_kegg_reactome)
            pickle.dump(self.__dict__, open(location_msigdb, 'wb'))


    def slurp_data(self, location):

        collected_pathways = []

        with open(location, 'r') as location_f:
            for line in location_f.readlines():
                line = line.strip()
                l = line.split('\t')
                pathway = Pathway(l[0], l[1], l[2:], self.entrez_id_to_symbol)
                collected_pathways.append(pathway)


        return collected_pathways



def main():

    msigdb = MsigDB()

    print msigdb.pathways_biocarta_kegg_reactome[0].entrez_gene_ids
    print msigdb.pathways_biocarta_kegg_reactome[0].gene_symbols
    print len(msigdb.pathways_biocarta_kegg_reactome)
    print len(msigdb.pathways_entrez)
    print msigdb.pathways_entrez[0].entrez_gene_ids
    print msigdb.pathways_entrez[0].gene_symbols
    print msigdb.entrez_id_to_symbol[127]


if __name__ == '__main__':
    main()

