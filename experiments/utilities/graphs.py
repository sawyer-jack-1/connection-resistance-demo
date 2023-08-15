# import numpy as np
import autograd.numpy as np

from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix
import networkx as nx
import matplotlib.pyplot as plt
import scipy

from utilities.rotation_matrices import Rxyz, Rx2D, Rx1D

class BaseGraph():
    def __init__(self, config):
        self.config = config.copy()

        self.edge_info = self.config['graph']['edge_info']
        self.edges = self.get_edges()
        self.numEdges = len(self.edges)

        self.numNodes = self.config['graph']['numNodes']
        self.nodes = self.get_nodes()

        self.weight_matrix = self.construct_weight_matrix()
        self.degree_matrix = self.construct_degree_matrix(self.weight_matrix)
        self.laplacian = self.construct_laplacian(self.degree_matrix, self.weight_matrix)
        self.pinv_laplacian = self.construct_pinv_laplacian(self.laplacian)

    def get_edges(self):
        edges = []
        for info in self.edge_info:
            edges.append((info[0], info[1]))
        return edges

    def get_nodes(self):
        return np.arange(0, self.numNodes)

    def get_laplacian(self):
        return self.laplacian

    def get_pinv_laplacian(self):
        return self.pinv_laplacian

    def get_edge_weights(self):
        edge_weights = []
        for info in self.edge_info:
            edge_weights.append(info[-1])
        return np.array(edge_weights)

    def construct_weight_matrix(self):
        weight_matrix = np.zeros((self.numNodes, self.numNodes))
        for edge_weight in self.config['graph']['edge_info']:
            i, j, w = edge_weight
            weight_matrix[i, j] = w
            weight_matrix[j, i] = w
        return weight_matrix

    def construct_degree_matrix(self, weight_matrix):
        return np.diag(np.sum(weight_matrix, axis=1))

    def construct_laplacian(self, degree_matrix, weight_matrix):
        return degree_matrix - weight_matrix

    def construct_pinv_laplacian(self, laplacian):
        return scipy.linalg.pinv(laplacian, rcond=10 ** (-10))


class ConnectionGraph(BaseGraph):
    def __init__(self, config):
        super().__init__(config)
        self.config = config.copy()
        self.edges = list(self.get_edges())
        self.nodes = list(self.get_nodes())
        self.numEdges = len(self.edges)
        self.numNodes = len(self.nodes)
        self.connectionDim = self.config['connectionGraph']['connectionDim']
        self.edgeRotations = self.config['connectionGraph']['edgeRotations']

        self.connection_incidence_matrix = self.construct_connection_incidence_matrix()
        self.connection_incidence_matrix_chung2014 = self.construct_connection_incidence_matrix_chung2014()
        self.connection_laplacian = self.construct_connection_laplacian()
        self.connection_laplacian_chung2014 = self.construct_connection_laplacian_chung2014()
        self.pinv_connection_laplacian = self.construct_pinv_connection_laplacian()
        self.pinv_connection_laplacian_chung2014 = self.construct_pinv_connection_laplacian_chung2014()

    #######################
    #### Get functions ####
    #######################
    def get_connectionDim(self):
        return self.connectionDim

    def get_numEdges(self):
        return self.numEdges

    def get_numNodes(self):
        return self.numNodes

    def get_connection_laplacian(self):
        return self.connection_laplacian

    def get_connection_laplacian_chung2014(self):
        return self.connection_laplacian_chung2014

    def get_connection_incidence_matrix(self):
        return self.connection_incidence_matrix

    def get_pinv_connection_laplacian(self):
        return self.pinv_connection_laplacian

    def get_pinv_connection_laplacian_chung2014(self):
        return self.pinv_connection_laplacian_chung2014

    def get_connection_incidence_matrix_chung2014(self):
        return self.connection_incidence_matrix_chung2014

    def get_connection_laplacian_chung2014(self):
        return self.connection_laplacian_chung2014

    def get_pinv_connection_laplacian_chung2014(self):
        return self.pinv_connection_laplacian_chung2014

    ############################
    #### Construct matrices ####
    #############################
    def construct_connection_incidence_matrix(self):
        connection_incidence_matrix = np.zeros((self.numEdges * self.connectionDim, self.numNodes * self.connectionDim))
        for edge_nr, edge_rot in enumerate(self.edgeRotations):
            e_start = edge_nr * self.connectionDim
            e_stop = (edge_nr + 1) * self.connectionDim
            i_start = self.edges[edge_nr][0] * self.connectionDim
            i_stop = (self.edges[edge_nr][0] + 1) * self.connectionDim
            j_start = self.edges[edge_nr][1] * self.connectionDim
            j_stop = (self.edges[edge_nr][1] + 1) * self.connectionDim

            connection_incidence_matrix[e_start:e_stop, i_start:i_stop] = np.diag(np.ones(self.connectionDim))
            if self.connectionDim == 1:
                connection_incidence_matrix[e_start:e_stop, j_start:j_stop] = - Rx1D(self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 2:
                connection_incidence_matrix[e_start:e_stop, j_start:j_stop] = - Rx2D(self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 3:
                connection_incidence_matrix[e_start:e_stop, j_start:j_stop] = - Rxyz(self.edgeRotations[edge_nr])
        return connection_incidence_matrix

    def construct_connection_laplacian(self):
        edge_weights = self.get_edge_weights()
        edge_weights_extended = np.diag(np.repeat(edge_weights, repeats=self.connectionDim))
        v = np.dot(edge_weights_extended, self.connection_incidence_matrix)
        return np.dot(self.connection_incidence_matrix.transpose(), v)

    def construct_pinv_connection_laplacian(self):
        return scipy.linalg.pinv(self.connection_laplacian, rcond=10 ** (-10))

    #######################################
    #### Construct matrices Chung 2014 ####
    #######################################

    def construct_connection_incidence_matrix_chung2014(self):
        connection_incidence_matrix_chung2014 = np.zeros(
            (self.numEdges * self.connectionDim, self.numNodes * self.connectionDim))
        for edge_nr, edge_rot in enumerate(self.edgeRotations):
            e_start = edge_nr * self.connectionDim
            e_stop = (edge_nr + 1) * self.connectionDim
            i_start = self.edges[edge_nr][0] * self.connectionDim
            i_stop = (self.edges[edge_nr][0] + 1) * self.connectionDim
            j_start = self.edges[edge_nr][1] * self.connectionDim
            j_stop = (self.edges[edge_nr][1] + 1) * self.connectionDim

            if self.connectionDim == 1:
                connection_incidence_matrix_chung2014[e_start:e_stop, i_start:i_stop] = Rx1D(
                    self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 2:
                connection_incidence_matrix_chung2014[e_start:e_stop, i_start:i_stop] = Rx2D(
                    self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 3:
                connection_incidence_matrix_chung2014[e_start:e_stop, i_start:i_stop] = Rxyz(self.edgeRotations[edge_nr])
            connection_incidence_matrix_chung2014[e_start:e_stop, j_start:j_stop] = - np.diag(np.ones(self.connectionDim))
        return connection_incidence_matrix_chung2014

    def construct_connection_laplacian_chung2014(self):
        edge_weights = self.get_edge_weights()
        edge_weights_extended = np.diag(np.repeat(edge_weights, repeats=self.connectionDim))
        v = np.dot(edge_weights_extended, self.connection_incidence_matrix_chung2014)
        return np.dot(self.connection_incidence_matrix_chung2014.transpose(), v)


    def construct_pinv_connection_laplacian_chung2014(self):
        return scipy.linalg.pinv(self.connection_laplacian_chung2014, rcond=10 ** (-10))
