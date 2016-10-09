#!/usr/bin/env python2.7



from load_networks import Network, NetworkBox
import glob, time
import numpy as np
import cPickle as pickle
import numpy as np
from load_msigdb import Pathway, MsigDB



msig_nets = MsigDB()

challenge_nets = NetworkBox()



#predictions1 = pickle.load(open('../../data/dream_networks/training/final/predictions1.pkl', 'rb'))

pred1_all = {}


#for net in predictions1:
#    p_net = {}
#    for cl in predictions1[net]:
#        p_net[cl] = predictions1[net][cl]
#    pred1_all[net] = p_net
#
#
#print 'done the 1st one'


subchallenge1_files = glob.glob('../../data/dream_networks/training/final/subchallenge1_*.pkl')

for sub1_f in subchallenge1_files:
    pred = pickle.load(open(sub1_f, 'rb'))
    cl = sub1_f.replace('../../data/dream_networks/training/final/subchallenge1_', '').replace('.pkl', '')
    for net in pred:
        pred1_all[net][cl] = pred[net]


print 'done 1st one'


print pred1_all[pred1_all.keys()[1]].keys()[:10]



thresholds = [.9, .8, .7, .6, .5, .4, .3, .2] + [10**(-x) for x in range(1, 32)]

for net in pred1_all.keys():
    indeces = range(len(pred1_all[net]))
    net_out_f = '../../data/dream_networks/model_application/subchallenge1/' + net
    net_out_f_handle = open(net_out_f, 'w')
    pathway_cntr = 1
    cl_done = []
    for thresh in thresholds:
        cl_count = 0
        for cl in pred1_all[net]:
            cl_count += 1
            if cl in cl_done:
                continue
            cl_len = [len(x.gene_symbols) for x in msig_nets.pathways_biocarta_kegg_reactome if x.name == cl][0]
            #print cl, cl_count, cl_len
            best = [i for i in indeces if pred1_all[net][cl][i] >= thresh]
            best = sorted(best, key = lambda i: pred1_all[net][cl][i])
            if len(best) >= cl_len:
                best = best[:cl_len]
                indeces = [x for x in indeces if x not in best]
                nodes = [str(challenge_nets.subchallenge1_networks[net].node_names[x]) for x in best]
                score = sum([pred1_all[net][cl][i] for i in best]) / len(best)
                net_out_f_handle.write(str(pathway_cntr) + '\t' + str(float(score)) + '\t' + '\t'.join(nodes) + '\n')
                cl_done.append(cl)
                pathway_cntr += 1
    net_out_f_handle.close()

print 'way done the 1st one'
print 'way done the 1st one for reals i think???'

print pred1_all.keys()

print pred1_all['1_ppi_anonym_v2.txt']['BIOCARTA_VEGF_PATHWAY'][pred1_all['1_ppi_anonym_v2.txt']['BIOCARTA_VEGF_PATHWAY'] > 0.01]








#predictions2 = pickle.load(open('../../data/dream_networks/training/final/predictions2.pkl', 'rb'))
#print challenge_nets.subchallenge2_networks['1_ppi_anonym_aligned_v2.txt'].node_names[:10]

all_genes = list(set(sum([challenge_nets.subchallenge2_networks[x].node_names for x in challenge_nets.subchallenge2_networks.keys()], [])))

all_pathways = [x.name for x in msig_nets.pathways_biocarta_kegg_reactome]

#print all_genes[:10]
#print challenge_nets.subchallenge2_networks['1_ppi_anonym_aligned_v2.txt'].node_names.index(0)

def get_gene_from_net(x, gene):
    try:
        return challenge_nets.subchallenge2_networks[x].node_names.index(gene)
    except ValueError:
        pass


positions = {}

for gene in all_genes:
    positions[gene] = {x : get_gene_from_net(x, gene) for x in challenge_nets.subchallenge2_networks.keys()}


part_pathways = predictions2['1_ppi_anonym_aligned_v2.txt'].keys()


#pred2_all = {cl : [] for cl in part_pathways}
pred2_all = {}

#for gene in all_genes:
#    p = positions[gene]
#    for cl in part_pathways:
#        pred2_all[cl].append(sum([float(predictions2[pk][cl][p[pk]]) for pk in p.keys() if p[pk] != None and p[pk] < predictions2[pk][cl].shape[0]]))


#print 'done'
#print predictions2[pk][cl]

subchallenge2_files = glob.glob('../../data/dream_networks/training/final/subchallenge2_*.pkl')

for sub2_f in subchallenge2_files:
    pred = pickle.load(open(sub2_f, 'rb'))
    cl = sub2_f.replace('../../data/dream_networks/training/final/subchallenge2_', '').replace('.pkl', '')
    pred2_all[cl] = [sum([float(pred[pk][positions[gene][pk]]) for pk in positions[gene].keys() if positions[gene][pk] != None and positions[gene][pk] < pred[pk].shape[0]]) for gene in all_genes]

print 'way done'


thresholds = [1.8, 1.5, 1.2, .9, .8, .7, .6, .5, .4, .3, .2] + [10**(-x) for x in range(1, 32)]

#print [max(pred2_all[cl]) for cl in pred2_all.keys()]

indeces = range(len(pred2_all[pred2_all.keys()[0]]))
net_out_f = '../../data/dream_networks/model_application/subchallenge2/sub2.txt'
net_out_f_handle = open(net_out_f, 'w')
pathway_cntr = 1
cl_done = []
for thresh in thresholds:
    cl_count = 0
    for cl in pred2_all:
        cl_count += 1
        if cl in cl_done:
            continue
        cl_len = [len(x.gene_symbols) for x in msig_nets.pathways_biocarta_kegg_reactome if x.name == cl][0]
        #print cl, cl_count, cl_len
        best = [i for i in indeces if pred2_all[cl][i] >= thresh]
        best = sorted(best, key = lambda i: pred2_all[cl][i])
        if len(best) >= cl_len:
            best = best[:cl_len]
            indeces = [x for x in indeces if x not in best]
            nodes = [str(all_genes[x]) for x in best]
            score = sum([pred2_all[cl][i] for i in best]) / len(best)
            net_out_f_handle.write(str(pathway_cntr) + '\t' + str(float(score)) + '\t' + '\t'.join(nodes) + '\n')
            cl_done.append(cl)
            pathway_cntr += 1
net_out_f_handle.close()


print 'all done'


