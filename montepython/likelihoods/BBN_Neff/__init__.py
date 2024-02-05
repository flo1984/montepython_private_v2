import os
import numpy as np
from montepython.likelihood_class import Likelihood
import scipy.constants as conts

class BBN_Neff(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)


        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        theo = cosmo.get_current_derived_parameters(['delta_N_eff_UV'])['delta_N_eff_UV'] + 3.044

        chi2 += ((theo - self.Neff) / self.error) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
