import os
from montepython.likelihood_class import Likelihood_prior


class OmegabBBN(Likelihood_prior):

    # initialisation of the class is done within the parent Likelihood_prior. For
    # this case, it does not differ, actually, from the __init__ method in
    # Likelihood class.
    def loglkl(self, cosmo, data):

        omegabvalue = cosmo.omega_b()
        loglkl = -0.5 * (omegabvalue - self.omegab) ** 2 / (self.sigma ** 2)

        return loglkl
