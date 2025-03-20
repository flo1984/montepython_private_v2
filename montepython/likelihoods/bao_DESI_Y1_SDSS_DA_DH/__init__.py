import os
import numpy as np
from montepython.likelihood_class import Likelihood
import montepython.io_mp as io_mp
import warnings


class bao_DESI_Y1_SDSS_DA_DH(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)


        # define array for values of z and data points
        self.z = np.array([], 'float64')
        self.dataDMrd = np.array([], 'float64')
        self.errorDMrd = np.array([], 'float64')
        self.dataDHrd = np.array([], 'float64')
        self.errorDHrd = np.array([], 'float64')
        self.corr = np.array([], 'float64')

        # read redshifts and data points
        with open(os.path.join(self.data_directory, self.file), 'r') as filein:
            for line in filein:
                if line.strip() and line.find('#') == -1:
                    # the first entry of the line is the identifier
                    this_line = line.split()
                    # insert into array if this id is not manually excluded
                    if not this_line[0] in self.exclude:
                        self.z = np.append(self.z, float(this_line[1]))
                        self.dataDMrd  = np.append(self.dataDMrd, float(this_line[2]))
                        self.errorDMrd = np.append(self.errorDMrd, float(this_line[3]))
                        self.dataDHrd  = np.append(self.dataDHrd, float(this_line[4]))
                        self.errorDHrd = np.append(self.errorDHrd, float(this_line[5]))
                        self.corr      = np.append(self.corr, float(this_line[6]))

        # number of data points
        self.num_points = np.shape(self.z)[0]

        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        # for each point, compute comoving distance DM and radial distance DH, and sound horizon at baryon drag rs_d, The quantities are already normalized to rs_d in the expression
        for i in range(self.num_points):
            
            DM = (1+self.z[i]) * cosmo.angular_distance(self.z[i]) /  cosmo.rs_drag()
            DH = 1. / cosmo.Hubble(self.z[i]) /  cosmo.rs_drag()
            
            # This just builds the covariance matrix
            Cov_mat = np.matmul(np.matmul(np.array([
                    [self.errorDMrd[i], 0],
                    [0, self.errorDHrd[i]]]),np.array([
                    [1, self.corr[i]],
                    [self.corr[i], 1]])),np.array([
                    [self.errorDMrd[i], 0],
                    [0, self.errorDHrd[i]]]))
                    
            invCov_mat = np.linalg.inv(Cov_mat)

            chi2 += np.dot(np.array([DM-self.dataDMrd[i],DH-self.dataDHrd[i]]),np.matmul(invCov_mat, np.array([DM-self.dataDMrd[i],DH-self.dataDHrd[i]])))
        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
        
