from utilities.rotation_matrices import Rxyz, Rx2D, Rx1D
from utilities.util import schur_comp
import numpy as np
from scipy import linalg

class ConnectionResistance():
    def __init__(self, config, connection_graph):
        self.connection_graph = connection_graph
        self.config = config

        self.connectionDim = self.config['connectionGraph']['connectionDim']
        self.numNodes = self.config['graph']['numNodes']

        self.edges = list(self.connection_graph.get_edges())

        self.pinv_laplacian = connection_graph.get_pinv_laplacian()

        self.connection_incidence_matrix = connection_graph.get_connection_incidence_matrix()
        self.connection_laplacian = connection_graph.get_connection_laplacian()
        self.pinv_connection_laplacian = connection_graph.get_pinv_connection_laplacian()

        self.connection_incidence_matrix_chung2014 = connection_graph.get_connection_incidence_matrix_chung2014()
        self.connection_laplacian_chung2014 = connection_graph.get_connection_laplacian_chung2014()
        self.pinv_connection_laplacian_chung2014 = connection_graph.get_pinv_connection_laplacian_chung2014()

    #######################
    #### Get functions ####
    #######################
    def get_edges(self):
        return self.edges

    ###################################
    #### Calc effective resistance ####
    ###################################
    def calc_effective_resistance(self):
        """Constructs the standard effective resistance matrix
        from the psudo inverse of the graph laplacian
        :return: Effective resistance matrix as n x n numpy matrix
        """
        n, n = self.pinv_laplacian.shape
        R = np.zeros((n, n))
        for i in range(0, n):
            for j in range(0, n):
                R[i, j] = self.pinv_laplacian[i, i] + self.pinv_laplacian[j, j] - 2 * self.pinv_laplacian[i, j]
        return R

    def calc_connection_resistance_ij(self, i, j):
        """Calculates the connection effective resistance between nodes i, j"""
        idx_i = self.connectionDim * i
        idx_j = self.connectionDim * j

        ones = np.ones(self.connectionDim)
        zeros = np.zeros(self.connectionDim)
        Xu = np.diag(np.concatenate((ones, zeros)))[:, 0:self.connectionDim]
        Xv = np.diag(np.concatenate((zeros, ones)))[self.connectionDim::, :].transpose()

        indices = np.array(
            [np.arange(idx_i, idx_i + self.connectionDim), np.arange(idx_j, idx_j + self.connectionDim)]).flatten()
        indices_comp = np.arange(0, self.numNodes * self.connectionDim)
        indices_comp = np.delete(indices_comp, indices, axis=0)

        sc = schur_comp(self.connection_laplacian, indices, indices_comp)
        c_ii = np.dot(Xu.transpose(), np.dot(sc, Xu))
        c_jj = np.dot(Xv.transpose(), np.dot(sc, Xv))

        c_ii_inv = np.linalg.inv(c_ii)
        c_jj_inv = np.linalg.inv(c_jj)

        r_ij = 1 / (2 * self.connectionDim) * np.trace(c_ii_inv + c_jj_inv)
        return r_ij

    def calc_connection_resistance(self):
        """Calculates the connection effective resistance matrix"""
        cr = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    cr[i, j] = self.calc_connection_resistance_ij(i, j)
        return cr

    def calc_connection_resistance_chung2014(self):
        """Calculates the connection effective resistance from Chung 20014"""
        psi = np.dot(np.dot(self.connection_incidence_matrix_chung2014, self.pinv_connection_laplacian_chung2014),
                     self.connection_incidence_matrix_chung2014.transpose())
        edges = self.connection_graph.get_edges()

        R = np.zeros((self.numNodes, self.numNodes))
        for k in range(0, len(edges)):
            i, j = edges[k]
            psi_ee = psi[k * self.connectionDim:(k + 1) * self.connectionDim,
                     k * self.connectionDim:(k + 1) * self.connectionDim]
            R[i, j] = np.linalg.norm(psi_ee, 2)
            R[j, i] = np.linalg.norm(psi_ee, 2)
        return R


