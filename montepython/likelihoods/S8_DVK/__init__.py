import os
import numpy as np
from montepython.likelihood_class import Likelihood
import scipy.constants as conts

class S8_DVK(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        self.need_cosmo_arguments(data, {'output': 'mPk'})

        if 'sigma8' not in data.get_mcmc_parameters(['derived']):
            raise io_mp.ConfigurationError('Error: S8 likelihood needs sigma8 as derived parameter')
        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        derived = cosmo.get_current_derived_parameters(data.get_mcmc_parameters(['derived']))
        for (name, value) in derived.items():
            data.mcmc_parameters[name]['current'] = value
        for name in derived:
            data.mcmc_parameters[name]['current'] /= data.mcmc_parameters[name]['scale']

        sigma8=data.mcmc_parameters['sigma8']['current']


        theo = sigma8*(cosmo.Omega_m()/0.3)**0.5
        print(theo)

        if theo-self.S8 > 0 :
            chi2 += ((theo - self.S8) / self.errorplus) ** 2
        else:
            chi2 += ((theo - self.S8) / self.errorminus) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
