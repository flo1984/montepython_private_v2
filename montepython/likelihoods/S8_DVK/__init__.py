import os
import numpy as np
from montepython.likelihood_class import Likelihood
import scipy.constants as conts

class S8_DVK(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)


        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        theo = cosmo.nl.sigma8[cosmo.nl.index_pk_m]*(cosmo.ba.Omega0_m/0.3)**0.5
        
        if theo-self.S8 > 0 :
            chi2 += ((theo - self.S8) / self.errorplus) ** 2
        else:
            chi2 += ((theo - self.S8) / self.errorminus) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
