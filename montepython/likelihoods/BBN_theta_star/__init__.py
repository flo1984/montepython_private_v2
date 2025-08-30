import os
import numpy as np
from montepython.likelihood_class import Likelihood
import scipy.constants as conts

class BBN_theta(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)


        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        omega_b_theo = cosmo.omega_b()
        theta_theo = (data.mcmc_parameters['M']['current'] *
             data.mcmc_parameters['M']['scale'])
        
        chi2 += ((omega_b_theo - self.omega_b) / self.error_omega_b) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
