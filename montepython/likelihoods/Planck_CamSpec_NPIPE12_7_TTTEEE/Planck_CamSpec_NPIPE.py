from montepython.likelihood_class import Likelihood
import numpy as np

def range_to_ells(use_range):
    """splits range string like '2-5 7 15-3000' into list of specific numbers"""
    if isinstance(use_range, str):
        ranges = []
        for ell_range in use_range.split():
            if '-' in ell_range:
                mn, mx = [int(x) for x in ell_range.split('-')]
                ranges.append(range(mn, mx + 1))
            else:
                ranges.append(int(ell_range))
        return np.concatenate(ranges)
    else:
        return use_range

class Planck_CamSpec_NPIPE(Likelihood):
    def __init__(self, path, data, command_line):
        Likelihood.__init__(self, path, data, command_line)
        ddir = self.data_directory
        spectra = np.loadtxt(ddir+self.cl_hat_file)
        if self.use_range is not None:
            used_ell = self.use_range
            if isinstance(used_ell, dict):
                print('Using range %s' % used_ell)
                used_ell = {key: range_to_ells(value) for key, value in used_ell.items()}
            else:
                print('CamSpec using range: %s' % used_ell)
                used_ell = range_to_ells(used_ell)
        else:
            used_ell = None
        data_vector = []
        nX = 0
        used_indices = []
        with open(ddir+self.data_ranges, "r", encoding="utf-8-sig") as f:
            lines = f.readlines()
            while not lines[-1].strip():
                lines = lines[:-1]
            self.Nspec = len(lines)
            lmin = np.zeros(self.Nspec, dtype=int)
            lmax = np.zeros(self.Nspec, dtype=int)
            self.cl_names = []
            self.ell_ranges = np.empty(self.Nspec, dtype=object)
            self.used_sizes = np.zeros(self.Nspec, dtype=int)
            for i, line in enumerate(lines):
                items = line.split()
                tp = items[0]
                self.cl_names.append(tp)
                lmin[i], lmax[i] = [int(x) for x in items[1:]]
                if lmax[i] and lmax[i] >= lmin[i]:
                    n = lmax[i] - lmin[i] + 1
                    data_vector.append(spectra[lmin[i]:lmax[i] + 1, i])
                    if tp in self.use_cl:
                        if used_ell is not None and (
                                not isinstance(used_ell, dict) or tp in used_ell):
                            if isinstance(used_ell, dict):
                                ells = used_ell[tp]
                            else:
                                ells = used_ell
                            self.ell_ranges[i] = np.array(
                                [L for L in range(lmin[i], lmax[i] + 1) if L in ells],
                                dtype=int)
                            used_indices.append(self.ell_ranges[i] + (nX - lmin[i]))
                        else:
                            used_indices.append(range(nX, nX + n))
                            self.ell_ranges[i] = range(lmin[i], lmax[i] + 1)
                        self.used_sizes[i] = len(self.ell_ranges[i])
                    else:
                        lmax[i] = -1
                    nX += n

        self.cl_used = np.array([name in self.use_cl for name in self.cl_names],
                                dtype=bool)
        covfile = ddir+self.covmat_fiducial
        with open(covfile, "rb") as cov_f:
            cov = np.fromfile(cov_f, dtype=[np.float32, np.float64]['64.bin' in covfile])
        assert (nX ** 2 == cov.shape[0])
        used_indices = np.concatenate(used_indices)
        self.data_vector = np.concatenate(data_vector)[used_indices]
        self.cov = cov.reshape(nX, nX)[np.ix_(used_indices, used_indices)].astype(
            np.float64)
        for name, mn, mx in zip(self.cl_names, lmin, lmax):
            if name in self.use_cl:
                print(name, mn, mx)
        print('Number of data points: %s' % self.cov.shape[0])
        self.lmax = lmax
        self.lmin = lmin
        max_l = np.max(self.lmax)
        self.ls = np.arange(max_l + 1)
        self.llp1 = self.ls * (self.ls + 1)
        self.covinv = np.linalg.inv(self.cov)

        self.need_cosmo_arguments(
            data, {'lensing': 'yes', 'output': 'tCl pCl lCl', 'l_max_scalars': self.lmax})

    def get_powerlaw_residuals(self, data_params):
        amp = np.empty(4)
        amp[0] = data_params['amp_100']['current']
        amp[1] = data_params['amp_143']['current']
        amp[2] = data_params['amp_217']['current']
        amp[3] = data_params['amp_143x217']['current']

        tilt = np.empty(4)
        tilt[0] = data_params['n_100']['current']
        tilt[1] = data_params['n_143']['current']
        tilt[2] = data_params['n_217']['current']
        tilt[3] = data_params['n_143x217']['current']

        powerlaw_pivot = 1500
        C_powerlaw = np.array(
            [amp[ii] * (self.ls / powerlaw_pivot) ** tilt[ii] for ii in range(4)])
        return C_powerlaw

    def get_cals(self, data_params):
        calPlanck = data_params['A_planck']['current'] ** 2
        calTE = data_params['calTE']['current']
        calEE = data_params['calEE']['current']
        return np.array([1, 1, 1, 1, calTE, calEE]) * calPlanck

    def chi_squared(self, CT, CTE, CEE, data_params):
        cals = self.get_cals(data_params)
        if np.any(self.cl_used[:4]):
            foregrounds = self.get_powerlaw_residuals(data_params)
        delta_vector = self.data_vector.copy()
        ix = 0
        for i, (cal, n) in enumerate(zip(cals, self.used_sizes)):
            if n > 0:
                if i <= 3:
                    # noinspection PyUnboundLocalVariable
                    delta_vector[ix:ix + n] -= (CT[self.ell_ranges[i]] +
                                                foregrounds[i][self.ell_ranges[i]]) / cal
                elif i == 4:
                    delta_vector[ix:ix + n] -= CTE[self.ell_ranges[i]] / cal
                elif i == 5:
                    delta_vector[ix:ix + n] -= CEE[self.ell_ranges[i]] / cal
                ix += n
        return np.dot(delta_vector, np.dot(self.covinv, delta_vector))

    def loglkl(self, cosmo, data):
        Cls = self.get_cl(cosmo)
        print(Cls)
        cltt, clte, clee, ll = Cls['tt'], Cls['te'], Cls['ee'], Cls['ell']
        fac = ll * (ll+1) / (2*np.pi)
        print(cltt,clte,clee,fac)
        cltt, clte, clee = cltt*fac, clte*fac, clee*fac
        data_params = data.mcmc_parameters
        lkl = -0.5*self.chi_squared(cltt, clte, clee, data_params)
        lkl = self.add_nuisance_prior(lkl, data)
        return lkl

    def add_nuisance_prior(self, lkl, data):
        # Recover the current value of the nuisance parameter.
        for nuisance in self.use_nuisance:
            nuisance_value = float(
                data.mcmc_parameters[nuisance]['current'] *
                data.mcmc_parameters[nuisance]['scale'])

            # add prior on nuisance parameters
            if hasattr(self, "%s_prior_center" % nuisance) and getattr(self, "%s_prior_std" % nuisance) > 0:
                # convenience variables
                prior_center = getattr(self, "%s_prior_center" % nuisance)
                prior_std = getattr(self, "%s_prior_std" % nuisance)
                lkl += -0.5*((nuisance_value-prior_center)/prior_std)**2
        return lkl
