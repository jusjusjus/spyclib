
import logging
import numpy as np
from .spaic2 import spaic2_betaspect as bs
from .spaic2 import spaic2_pot2density as p2d

logging.basicConfig(level=logging.INFO)
EPS = 1.0**-6


class Spaic2Solver:

    logger = logging.getLogger(name='Spaic2Solver')

    default_woodsaxon_params = np.array([
            0.84, 0.39, 0.78, 0.80, 0.91, 0.20, 0.34, 0.77, 0.28, 0.55
    ])

    default_quantum_numbers = np.array([
            [2, 1]
    ], dtype=np.int32)
    
    default_cluster_radius = np.float64(11.)# rclus
    default_cluster_steepness = np.float64(0.9) # drclus
 
    def __init__(self):
        self.logger.debug('__init__')
        # Define flags prior to calling self.init()
        self.cluster_steepness_flag = False
        self.cluster_radius_flag = False
        self.woodsaxon_params_flag = False
        self.init()
        self.spectra_computed = False
        

    def init(self):
        """Initialize the default potential and qms

        @property'es initial point to undefined states, thus
        necessitating an alternative initialization.
        """
        self.cluster_radius = self.default_cluster_radius # setter
        self.cluster_steepness = self.default_cluster_steepness # setter
        self.woodsaxon_params  = self.default_woodsaxon_params # setter
        qns = self.default_quantum_numbers
        # Here p2d is accessed, not bs.
        prs = np.array(p2d.pararho, dtype=np.float64)
        is_valid = len(qns.shape) == 2 and qns.shape[1] == 2
        assert is_valid, "QNs example: {} (received {})".format(self.default_quantum_numbers, qns)
        bs.set_density_and_qns(prs, qns)


    def plot(self, **kwargs):
        if not self.spectra_computed:
            self.compute_spectra()

        show = kwargs.get('show', True)
        fs = kwargs.get('fontsize', 14)
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 10))
        ax = plt.subplot(311)
        plt.plot(p2d.rr, self.woodsaxon_potential, label='Wood-Saxon Pot.')
        plt.plot(bs.rr, self.density_potential, label='e- Density Pot.')
        plt.legend(loc=0)
        plt.grid()
        ax = plt.subplot(312)
        for q, qn in enumerate(self.quantum_numbers):
            plt.plot(self.frequencies, self.spectra[:, q], label=str(qn))
        plt.grid()
        plt.subplot(313, sharex=ax)
        for q, qn in enumerate(self.quantum_numbers):
            plt.plot(self.frequencies, self.beta_spectra[:, q], label=str(qn))
        plt.grid()
        plt.xlim(self.frequencies.min(), self.frequencies.max())

        if show: plt.show()

    # Common interface
    def __del__(self):
        self.logger.debug('__del__')
        p2d.free_all_memory()
        bs.free_all_memory()

    @property
    def cluster_radius(self):
        """Return cluster radius used for generating the charge distribution.
        """
        return p2d.rclus

    @cluster_radius.setter
    def cluster_radius(self, params):
        """Set the cluster radius used for generating the charge distribution.
        """
        params = np.float64(params)
        assert isinstance(params,np.float64), "received {}".format(params)
        # set variable in p2d
        #p2d.rclus = params
        p2d.rclus = np.array(params, dtype=np.float64)
        self.cluster_radius_flag = True
        # Determine the density parameters and transmit them to betaSpect
        if self.cluster_steepness_flag and self.woodsaxon_params_flag:
            # "If condition" only important for __init__ where cluster_steepness is not set.
            p2d.compute_density()
            self.update_density_params(p2d.pararho)
    
    @property
    def cluster_steepness(self):
        """Return cluster steepness used for generating the charge distribution.
        """
        return p2d.drclus
    
    @cluster_steepness.setter
    def cluster_steepness(self, params):
        """Set the cluster steepness used for generating the charge distribution.
        """
        params = np.float64(params)
        assert isinstance(params,np.float64), "received {}".format(params)
        # set variable in p2d
        p2d.drclus = params
        self.cluster_steepness_flag = True
        # Determine the density parameters and transmit them to betaSpect
        if self.cluster_radius_flag and self.woodsaxon_params_flag:
            p2d.compute_density()
            self.update_density_params(p2d.pararho)
    
    @property
    def density_params(self):
        """Return c coefficients describing the charge distribtion.
        """
        return bs.para

    @density_params.setter
    def density_params(self, params):
        params = np.asarray(params, dtype=np.float64)
        assert len(params.shape) == 1, "received {}".format(params)
        self.update_density_params(params)
    

    @property
    def woodsaxon_params(self):
        """Return B coefficients defining the bumped Wood-Saxon potential.
        """
        return 0.5 * (p2d.parapot + 1.0)

    @woodsaxon_params.setter
    def woodsaxon_params(self, params):
        """Set B coefficients defining the bumped Wood-Saxon potential.
        """
        self.logger.debug('woodsaxon_params.setter')
        params = np.asarray(params, dtype=np.float64)
        assert len(params.shape) == 1, "received {}".format(params)
        assert all(params >= 0.0) and all(params <= 1.0), "received {}".format(params)
        p2d.set_woodsaxon_params(params)
        self.woodsaxon_params_flag = True
        # Compute new p2d.pararho (!= density_params which is bs.para)
        if self.cluster_radius_flag and self.cluster_steepness_flag:
            # "If condition" only important for __init__ where cluster 
            # parameters might not yet be set.
            self.logger.debug('compute_density')
            # p2d.diagnostics()
            p2d.compute_density()
            # Transmit pararho from pot2density to betaSpect
            self.update_density_params(p2d.pararho)

    @property
    def woodsaxon_potential(self):
        """Return the bumped Wood-Saxon potential.
        """
        return p2d.woodsaxon_potential

    # Interface to spaic2_betaspect
    @property
    def beta_spectra(self):
        assert self.spectra_computed, "Explicitely call compute_spectra first."
        return bs.beta

    @property
    def spectra(self):
        assert self.spectra_computed, "Explicitely call compute_spectra first."
        return bs.spec

    @property
    def frequencies(self):
        return bs.freq

    @property
    def density_potential(self):
        assert self.spectra_computed, "Explicitely call compute_spectra first."
        return bs.pot
    
    @property
    def radius(self):
        assert self.spectra_computed, "Explicitely call compute_spectra first."
        return bs.rr
    
    @property
    def bpsi(self):
        assert self.spectra_computed, "Explicitely call compute_spectra first."
        return bs.bpsi+bs.eb

    @property
    def quantum_numbers(self):
        return bs.quantum_numbers

    @quantum_numbers.setter
    def quantum_numbers(self, qns):
        qns = np.array(qns, dtype=np.int32)
        prs = np.array(self.density_params, dtype=np.float64)
        is_valid = len(qns.shape) == 2 and qns.shape[1] == 2
        assert is_valid, "QNs example: {} (received {})".format(self.default_quantum_numbers, qns)
        bs.set_density_and_qns(prs, qns)

    def update_density_params(self, dparams):
        self.logger.debug('update_density_params')
        prs = np.array(dparams)
        self.spectra_computed = False
        # Only important for __init__ where quantum_numbers are None.
        if self.quantum_numbers is not None:
            qns = np.array(self.quantum_numbers)
            bs.set_density_and_qns(prs, qns)

    def compute_spectra(self):
        bs.compute_beta()
        self.spectra_computed = True
