#!/usr/bin/env python2.7



import os.path
import cPickle as pickle
from glob import glob
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd




class Network:

    def __init__(self, network_location, **kwargs):
        """
        Networks are loaded here as Pandas objects, and transformed into sparse matrices that can be
        unpacked into full matrices. For unidirectional networks the matrix is non-symmetric, for bi-
        directional ones it's symmetric.

        Note on network names: these are guaranteed ot be in the correct order, w.r.t. the full rows
        and columns of the full matrix.
        """


        def sparse_transform(df):
            """
            This creates sparse matrices, in a surprisingly fast fashion. There's this backwards
            way of handling missing values here, where they get substituted with 0's if you unroll
            the sparse data into dense data with todense(), but apparently there is nothing that can
            be done with it and this behavior is deeply ingrained into NumPy. It should be ok here,
            as the values in all the networks here are never 0.
            """

            df_row_u = np.sort(df['node1'].unique()).tolist()
            df_col_u = np.sort(df['node2'].unique()).tolist()

            data = df['weight'].tolist()

            row = df['node1'].astype('category', categories = df_row_u).cat.codes
            col = df['node2'].astype('category', categories = df_col_u).cat.codes
            sparse_matrix = csr_matrix((data, (row, col)), shape = (len(df_row_u), len(df_col_u)))

            return sparse_matrix



        def symmetricize(x, y):
            """
            Transforms the data in dataframes so that they will produce symmetric sparse dataframes at a later step.
            """

            if x in ['3_signal_anonym_directed_v3.txt', '3_signal_anonym_aligned_directed_v3.txt']:
                pass

            else:
                y_pt2 = y[['node2', 'node1', 'weight']]
                y_pt2.columns = ['node1', 'node2', 'weight']
                y = pd.concat([y, y_pt2])

            return y



        if 'network_name' in kwargs.keys():
            network_name = kwargs['network_name']
        else:
            network_name = os.path.basename(network_location)


        if 'sep' in kwargs.keys():
            sep = kwargs['sep']
        else:
            sep = '\t'



        network_weights = pd.read_csv(network_location, sep = sep, index_col = False, header = None, names = ['node1', 'node2', 'weight'])
        network_weights = symmetricize(network_name, network_weights)
        self.network_weights_original = network_weights
        self.network_weights = sparse_transform(network_weights)
        self.node_names = np.sort(network_weights['node1'].append(network_weights['node2']).unique()).tolist()



    def transform_network(self):
        """
        'weight' values are taken directly from the challenge, and they are either confidence scores, or correlations.
        As such they can be used directly as some sort of force constants in the annealing. This function is a place-
        holder to change them to a different form before proceeding further, if the core values end up not working too
        well.
        """
        pass





class NetworkBox:
    """
    This class loads all networks from subchallenge 1 and 2. It then transforms the node-pair
    data to a distance space (long-to-wide transform), in various forms (best option needs to
    be picked).
    """

    def __init__(self):


        location_networks = '../../data/dream_networks/networks.pkl'

        if os.path.isfile(location_networks):
            temp_dict = pickle.load(open(location_networks, 'rb'))
            self.__dict__.update(temp_dict)

        else:

            subchallenge1_network_files = glob('../../data/dream_networks/subchallenge1/*_*.txt')
            subchallenge2_network_files = glob('../../data/dream_networks/subchallenge2/*_*.txt')


            self.subchallenge1_networks = {os.path.basename(x) : Network(x) for x in subchallenge1_network_files}
            self.subchallenge2_networks = {os.path.basename(x) : Network(x) for x in subchallenge2_network_files}

            pickle.dump(self.__dict__, open(location_networks, 'wb'))





def main():

    networks = NetworkBox()
    #print 'done'
    #print networks.subchallenge1_networks.keys()
    #print networks.subchallenge1_networks['1_ppi_anonym_v2.txt'].network_weights
    #print networks.subchallenge1_networks['1_ppi_anonym_v2.txt'].network_weights[17380, 17393]
    #print networks.subchallenge1_networks['1_ppi_anonym_v2.txt'].network_weights[17393, 17380]



if __name__ == '__main__':

    main()


