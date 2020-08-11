import os
import numpy as np
from montepython.likelihood_class import Likelihood
import scipy.constants as conts

class prior_omega_b(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)


        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        theo = cosmo.omega_b()

        chi2 += ((theo - self.omega_b) / self.error) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
