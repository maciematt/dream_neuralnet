#!/usr/bin/env python2.7



from load_msigdb import Pathway, MsigDB
from load_networks import Network, NetworkBox
from collections import defaultdict, Counter
from subprocess import call
from operator import itemgetter
import cPickle as pickle
import pandas as pd
import numpy as np
import time, gzip, os, random, itertools
import theano


from keras.wrappers.scikit_learn import KerasClassifier
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.utils.np_utils import to_categorical
from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import cross_val_score




msig_db = MsigDB()


gene_path_bkr = defaultdict(list)

for pathway in msig_db.pathways_biocarta_kegg_reactome:
    #if pathway.name in collated_test_scorelist_top_qval['Name'].tolist(): # if using only GWAS-significant pathways
    #if pathway.name in collated_test_scorelist['Name'].unique().tolist():
    for gene in pathway.gene_symbols:
        gene_path_bkr[gene].append(pathway.name)


genes = gene_path_bkr.keys()


#gene_string_significant = '../../data/gene_string_significant.txt' # if using only GWAS-significant pathways
#gene_symbols_significant = '../../data/gene_symbols_significant.txt' # if using only GWAS-significant pathways
gene_string_significant = '../../data/gene_string_allmsigdb.txt'
gene_symbols_significant = '../../data/gene_symbols_allmsigdb.txt'


if not os.path.isfile(gene_string_significant):

    with open(gene_symbols_significant, 'w') as gene_writer:
        gene_writer.write('%0D'.join(genes))

    call(['./map_gene_ids_to_string.sh', gene_symbols_significant, gene_string_significant])



gene_string_df = pd.read_csv(gene_string_significant, sep = '\t')
gene_string_path_bkr = {}
for gene in gene_path_bkr:
    string_id = gene_string_df[gene_string_df['preferredName'] == gene]['stringId'].tolist()
    if len(string_id) > 0:
        gene_string_path_bkr[string_id[0]] = gene_path_bkr[gene]



significant_string_genes = gene_string_path_bkr.keys()


full_string_file = '../../data/protein.links.v10.txt.gz'
#significant_string_file ='../../data/significant_part_string.txt' # if using only GWAS-significant pathways
significant_string_file ='../../data/allmsigdb_part_string.txt'


if not os.path.isfile(significant_string_file):

    chunksize = 1 * 10**8

    i = 0
    for chunk in pd.read_csv(full_string_file, sep = ' ', chunksize = chunksize):
        i += 1
        start = time.time()
        print 'Working on chunk ', i, '...'
        chunk_trimmed = chunk[chunk['protein1'].isin(significant_string_genes) & chunk['protein2'].isin(significant_string_genes)]
        chunk_trimmed.to_csv(significant_string_file, index = False, header = False, mode = 'a', sep = ' ')




string_network_file = '../../data/dream_networks/training/string/string_ready.pkl'



if os.path.isfile(string_network_file):

    string_network = pickle.load(open(string_network_file, 'rb'))

else:

    string_network = Network(significant_string_file, sep = ' ')


    string_network.network_mat = string_network.network_weights.todense()
    network_mat_mean = string_network.network_mat[string_network.network_mat != 0].mean()
    network_mat_sd = string_network.network_mat[string_network.network_mat != 0].std()
    network_mat_min = ((string_network.network_mat[string_network.network_mat != 0] - network_mat_mean) / network_mat_sd).min()
    network_mat_max = ((string_network.network_mat[string_network.network_mat != 0] - network_mat_mean) / network_mat_sd).max()


    # NEXT, for every row extract ONLY non-zero values and normalize them. Then create for every row (/sample/gene) a sample of the distribution of normalized values - as big a sample as the number of input dimensions I want, like say 100 samples so that I have 100 dims.


    # Another thing I could do is calculate the fingerprint of each node... this could be pretty nice but let's see if it's necessary?


    def draw_from_quantile(row):
        r = (row[row != 0 ] - network_mat_mean) / network_mat_sd
        r_out = np.histogram(r, range = (network_mat_min, network_mat_max), bins = 1000, density = True)
        return r_out


    string_network.network_mat_resampl = np.apply_along_axis(draw_from_quantile, 1, string_network.network_mat)[:, 0]

    print string_network.network_mat_resampl[0].shape


    string_network.network_mat_full = np.empty((0, 1000))
    string_network.node_names_full = []
    string_network.class_full = []


    for i in range(len(string_network.node_names)):
        if i % 100 == 0:
            print 'Now on i = ', i
        gene = string_network.node_names[i]
        classes = list(set(gene_string_path_bkr[gene]))
        string_network.class_full += classes
        string_network.node_names_full += [gene] * len(classes)
        row = string_network.network_mat_resampl[i]
        row = row[None, :]
        row_repl = np.tile(row, (len(classes), 1))
        string_network.network_mat_full = np.append(string_network.network_mat_full, row_repl, axis = 0)



    pickle.dump(string_network, open(string_network_file, 'wb'))





def upsample(positive_indeces):
    """
    """

    rand_indeces = list(set(range(len(string_network.class_full))) - set(positive_indeces))

    resample_minor = len(rand_indeces) / len(positive_indeces)

    rest = len(rand_indeces) - resample_minor * len(positive_indeces)
    rest_resampled = random.sample(positive_indeces, rest)

    positive_array = np.tile(string_network.network_mat_full[positive_indeces, :], reps = (resample_minor, 1))
    positive_array = np.append(positive_array, string_network.network_mat_full[rest_resampled, :], axis = 0)

    full_array = np.append(positive_array, string_network.network_mat_full[rand_indeces, :], axis = 0)

    y_values = [1] * positive_array.shape[0] + [0] * positive_array.shape[0]


    return full_array, y_values



def downsample(positive_indeces):
    """
    """

    rand_indeces = list(set(range(len(string_network.class_full))) - set(positive_indeces))
    rest_resampled = random.sample(rand_indeces, len(positive_indeces))

    full_array = np.append(string_network.network_mat_full[positive_indeces, :], string_network.network_mat_full[rest_resampled, :], axis = 0)

    y_values = [1] * len(positive_indeces) + [0] * len(positive_indeces)

    return full_array, y_values





#def create_model():
#    """
#    Not bad - avg acc 0.71
#    """
#
#    model_bin = Sequential()
#
#    model_bin.add(Dense(output_dim = 1000, input_dim = 1000))
#    model_bin.add(Activation("relu"))
#    model_bin.add(Dense(output_dim = 1250, input_dim = 1000))
#    model_bin.add(Activation("relu"))
#    model_bin.add(Dense(output_dim = 750, input_dim = 1250))
#    model_bin.add(Activation("relu"))
#    #model_bin.add(Dense(output_dim = 500, input_dim = 750))
#    #model_bin.add(Activation("relu"))
#    #model_bin.add(Dense(output_dim = 750, input_dim = 500))
#    #model_bin.add(Activation("relu"))
#    model_bin.add(Dense(output_dim = 1))
#    model_bin.add(Activation("sigmoid"))
#
#    model_bin.compile(loss = 'binary_crossentropy', optimizer = 'adam', metrics = ['accuracy'])
#
#    return model_bin



class ResultObj:
    def __init__(self, results, results_mean, model_params):
        self.results = results
        self.results_mean = results_mean
        self.model_params = model_params




def run_model_variant(var_num):

    model_params = dict(zip(['layer_1', 'layer_2', 'layer_3', 'layer_4', 'input_dropout', 'hidden_dropout'], pos_permutations[var_num]))


    def create_model():
        """
        """

        model_bin = Sequential()

        model_bin.add(Dense(output_dim = model_params['layer_1'], input_dim = 1000))
        model_bin.add(Activation("relu"))

        if model_params['input_dropout']:
            model_bin.add(Dropout(0.2))

        if model_params['layer_2']:
            model_bin.add(Dense(output_dim = model_params['layer_2'], input_dim = model_params['layer_1']))
            model_bin.add(Activation("relu"))

        if model_params['hidden_dropout']:
            model_bin.add(Dropout(0.5))

        if model_params['layer_3']:
            model_bin.add(Dense(output_dim = model_params['layer_3'], input_dim = model_params['layer_2']))
            model_bin.add(Activation("relu"))

        if model_params['layer_4']:
            model_bin.add(Dense(output_dim = model_params['layer_4'], input_dim = model_params['layer_3']))
            model_bin.add(Activation("relu"))

        model_bin.add(Dense(output_dim = 1))
        model_bin.add(Activation("sigmoid"))

        model_bin.compile(loss = 'binary_crossentropy', optimizer = 'adam', metrics = ['accuracy'])

        return model_bin



    model_bin = KerasClassifier(build_fn = create_model, nb_epoch = 15, batch_size = 16)

    seed = 123
    kfold = StratifiedKFold(y = y_values, n_folds = 5, shuffle = True, random_state = seed)


    results = cross_val_score(model_bin, full_array, y_values, cv = kfold)
    result_object = ResultObj(results, results.mean(), model_params)

    return result_object


#model_bin.fit(full_array, y_values, nb_epoch = 40, batch_size = 16)




result_list = []

def log_result(result):
    r = ResultObj(*result)
    result_list.append(r)
    print 'So far done', len(result_list), ' out of', len(pos_permutations)



positive_indeces = [i for i, e in enumerate(string_network.class_full) if e == string_network.class_full[0]]
full_array, y_values = downsample(positive_indeces)


#layer_1 = [250, 500, 1000, 2000, 4000, 8000, 16000]
#layer_2 = [False, 250, 500, 1000, 2000, 4000, 8000, 16000]
#layer_3 = [False, 250, 500, 1000, 2000, 4000, 8000, 16000]
#layer_4 = [False, 250, 500, 1000, 2000, 4000, 8000, 16000]
#input_dropout = [False, True]
#hidden_dropout = [False, True]

layer_1 = [1000, 2000, 4000]
layer_2 = [False, 1000, 2000, 4000]
layer_3 = [False, 1000, 2000, 4000]
layer_4 = [False, 1000, 2000, 4000]
input_dropout = [False]
hidden_dropout = [False]


combined_possibilities = [layer_1, layer_2, layer_3, layer_4, input_dropout, hidden_dropout]
pos_permutations = list(itertools.product(*combined_possibilities))


result_list = []
for var_num in range(len(pos_permutations)):
    print '\nworking on', var_num + 1, ' out of', len(pos_permutations)
    result = run_model_variant(var_num)
    result_list.append(result)


pickle.dump(result_list, open('../../data/dream_networks/training/string/cv_test_result_list.pkl', 'wb'))

best_index = max(enumerate([r.results_mean for r in result_list]), key = itemgetter(1))[0]


layer_1_2 = [250, 500, 1000]
layer_2_2 = [250, 500, 1000, 2000]
layer_3_2 = [250, 500, 1000, 2000]
layer_4_2 = [False]
combined_possibilities_2 = [layer_1_2, layer_2_2, layer_3_2, input_dropout, hidden_dropout]
pos_permutations_2 = list(itertools.product(*combined_possibilities_2))


result_list_2 = []
for var_num in range(len(pos_permutations_2)):
    print '\nworking on', var_num + 1, ' out of', len(pos_permutations_2)
    result = run_model_variant(var_num)
    result_list_2.append(result)

best_index_2 = max(enumerate([r.results_mean for r in result_list_2]), key = itemgetter(1))[0]



layer_1_3 = [500, 750, 1000]
layer_2_3 = [500, 750, 1000, 2000]
layer_3_3 = [500, 750, 1000]
layer_4_3 = [False]
combined_possibilities_3 = [layer_1_3, layer_2_3, layer_3_3, input_dropout, hidden_dropout]
pos_permutations_3 = list(itertools.product(*combined_possibilities_3))


result_list_3 = []
for var_num in range(len(pos_permutations_3)):
    print '\nworking on', var_num + 1, ' out of', len(pos_permutations_3)
    result = run_model_variant(var_num)
    result_list_3.append(result)


best_index_3 = max(enumerate([r.results_mean for r in result_list_3]), key = itemgetter(1))[0]


combined_possibilities_4 = [[[1000, 2000, 2000], [500, 2000, 500], [500, 500, 1000]], [False], [False, True], [False, True]]
pos_permutations_4 = list(itertools.product(*combined_possibilities_4))
pos_permutations_4 = [[e[0][0], e[0][1], e[0][2], e[1], e[2], e[3]] for e in pos_permutations_4]

result_list_4 = []
for var_num in range(len(pos_permutations_4)):
    print '\nworking on', var_num + 1, ' out of', len(pos_permutations_4)
    result = run_model_variant(var_num)
    result_list_4.append(result)


best_index_4 = max(enumerate([r.results_mean for r in result_list_4]), key = itemgetter(1))[0]
print result_list_4[best_index_4].results_mean


final_parameters = pos_permutations_4[best_index_4]


pickle.dump([result_list, result_list_2, result_list_3, result_list_4, final_parameters], open('../../data/dream_networks/training/string/cv_test_result_list.pkl', 'wb'))

#result_list, result_list_2, result_list_3, result_list_4, final_parameters = pickle.load(open('../../data/dream_networks/training/string/cv_test_result_list.pkl', 'rb'))




def create_final_model():
    """
    """

    model_bin = Sequential()

    model_bin.add(Dense(output_dim = final_parameters[0], input_dim = 1000))
    model_bin.add(Activation("relu"))

    model_bin.add(Dropout(0.2))

    model_bin.add(Dense(output_dim = final_parameters[1], input_dim = final_parameters[0]))
    model_bin.add(Activation("relu"))

    model_bin.add(Dropout(0.5))

    model_bin.add(Dense(output_dim = final_parameters[2], input_dim = final_parameters[1]))
    model_bin.add(Activation("relu"))


    model_bin.add(Dense(output_dim = 1))
    model_bin.add(Activation("sigmoid"))

    model_bin.compile(loss = 'binary_crossentropy', optimizer = 'adam', metrics = ['accuracy'])

    return model_bin



class PrunedNet:
    def __init__(self):
        pass


def draw_from_quantile(row):
    r = (row[row != 0 ] - network_mat_mean) / network_mat_sd
    r_out = np.histogram(r, range = (network_mat_min, network_mat_max), bins = 1000, density = True)
    return r_out



dream_networks = NetworkBox()


dream_networks_for_calcs = {}
dream_networks_for_calcs['subchallenge1'] = {}
dream_networks_for_calcs['subchallenge2'] = {}
predictions1 = {}
predictions2 = {}


for net in dream_networks.subchallenge1_networks:
    dream_networks.subchallenge1_networks[net].network_mat = dream_networks.subchallenge1_networks[net].network_weights.todense()
    network_mat_mean = dream_networks.subchallenge1_networks[net].network_mat[dream_networks.subchallenge1_networks[net].network_mat != 0].mean()
    network_mat_sd = dream_networks.subchallenge1_networks[net].network_mat[dream_networks.subchallenge1_networks[net].network_mat != 0].std()
    network_mat_min = ((dream_networks.subchallenge1_networks[net].network_mat[dream_networks.subchallenge1_networks[net].network_mat != 0] - network_mat_mean) / network_mat_sd).min()
    network_mat_max = ((dream_networks.subchallenge1_networks[net].network_mat[dream_networks.subchallenge1_networks[net].network_mat != 0] - network_mat_mean) / network_mat_sd).max()

    temp = np.apply_along_axis(draw_from_quantile, 1, dream_networks.subchallenge1_networks[net].network_mat)[:, 0]
    x_test = np.array([temp[x, ] for x in range(temp.shape[0])])

    dream_networks_for_calcs['subchallenge1'][net] = PrunedNet()

    dream_networks_for_calcs['subchallenge1'][net].network_mat_resampl = x_test
    dream_networks_for_calcs['subchallenge1'][net].node_names = dream_networks.subchallenge1_networks[net].node_names
    predictions1[net] = {}


del dream_networks.subchallenge1_networks


for net in dream_networks.subchallenge2_networks:
    dream_networks.subchallenge2_networks[net].network_mat = dream_networks.subchallenge2_networks[net].network_weights.todense()
    network_mat_mean = dream_networks.subchallenge2_networks[net].network_mat[dream_networks.subchallenge2_networks[net].network_mat != 0].mean()
    network_mat_sd = dream_networks.subchallenge2_networks[net].network_mat[dream_networks.subchallenge2_networks[net].network_mat != 0].std()
    network_mat_min = ((dream_networks.subchallenge2_networks[net].network_mat[dream_networks.subchallenge2_networks[net].network_mat != 0] - network_mat_mean) / network_mat_sd).min()
    network_mat_max = ((dream_networks.subchallenge2_networks[net].network_mat[dream_networks.subchallenge2_networks[net].network_mat != 0] - network_mat_mean) / network_mat_sd).max()

    temp = np.apply_along_axis(draw_from_quantile, 1, dream_networks.subchallenge2_networks[net].network_mat)[:, 0]

    x_test = np.array([temp[x, ] for x in range(temp.shape[0])])

    dream_networks_for_calcs['subchallenge2'][net] = PrunedNet()

    dream_networks_for_calcs['subchallenge2'][net].network_mat_resampl = x_test
    dream_networks_for_calcs['subchallenge2'][net].node_names = dream_networks.subchallenge2_networks[net].node_names
    predictions2[net] = {}


del dream_networks.subchallenge2_networks



unique_classes = list(set(string_network.class_full))

final_models = {}


mod_num = 0

start = time.time()
for cl in unique_classes:

    mod_num += 1

    #if mod_num < 110:
    #    continue

    print 'Working on model', mod_num, ' out of', len(unique_classes)


    positive_indeces = [i for i, e in enumerate(string_network.class_full) if e == cl]
    full_array, y_values = upsample(positive_indeces)

    model_bin = create_final_model()

    model_bin.fit(full_array, y_values, nb_epoch = 3, batch_size = 2048, verbose = 0)
    #final_models[cl] = model_bin

    print 'Done!', (time.time() - start) / 60, ' min elapsed so far...'

    pred1 = {}
    for net in dream_networks_for_calcs['subchallenge1']:
        pred1[net] = model_bin.predict(dream_networks_for_calcs['subchallenge1'][net].network_mat_resampl)
        #predictions1[net][cl] = pred

    pred2 = {}
    for net in dream_networks_for_calcs['subchallenge2']:
        pred2[net] = model_bin.predict(dream_networks_for_calcs['subchallenge2'][net].network_mat_resampl)
        #predictions2[net][cl] = pred



    pickle.dump(pred1, open('../../data/dream_networks/training/final/subchallenge1_' + cl + '.pkl', 'wb'))
    pickle.dump(pred2, open('../../data/dream_networks/training/final/subchallenge2_' + cl + '.pkl', 'wb'))




#pickle.dump(final_models, open('../../data/dream_networks/training/string/models_fit.pkl', 'wb'))


#cl = 'BIOCARTA_VEGF_PATHWAY'
#positive_indeces = [i for i, e in enumerate(string_network.class_full) if e == cl]
#full_array, y_values = upsample(positive_indeces)
#test = final_models[cl].predict(full_array)
print 'done'
#print final_models

