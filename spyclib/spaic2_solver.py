
import logging
import numpy as np
from .spaic2 import spaic2_betaspect as bs
from .spaic2 import spaic2_d2pot as d2p

logging.basicConfig(level=logging.INFO)
EPS = 1.0**-6


class Spaic2Solver:

    logger = logging.getLogger(name='Spaic2Solver')

    default_woodsaxon_params = np.array([
            0.84, 0.39, 0.78, 0.80, 0.91, 0.20, 0.34, 0.77, 0.28,
            0.55, 0.48, 0.63, 0.36, 0.51, 0.95, 0.92, 0.64, 0.72, 0.14, 0.61,
            0.02, 0.24, 0.14, 0.80, 0.16, 0.40, 0.13, 0.11, 1.00, 0.22
    ])

    default_quantum_numbers = np.array([
            [0, 1],
            [1, 1]
    ], dtype=np.int32)

    def __init__(self):
        self.woodsaxon_params   = self.default_woodsaxon_params # setter
        self.quantum_numbers  = self.default_quantum_numbers # setter
        self.spectra_computed = False

    def plot(self, **kwargs):
        if not self.spectra_computed:
            self.compute_spectra()

        show = kwargs.get('show', True)
        fs = kwargs.get('fontsize', 14)
        import matplotlib.pyplot as plt

        ax = plt.subplot(211)
        for q, qn in enumerate(self.quantum_numbers):
            plt.plot(self.frequencies, self.spectra[:, q], label=str(qn))
        plt.grid()
        plt.subplot(212, sharex=ax)
        for q, qn in enumerate(self.quantum_numbers):
            plt.plot(self.frequencies, self.beta_spectra[:, q], label=str(qn))
        plt.grid()
        plt.xlim(self.frequencies.min(), self.frequencies.max())

        if show: plt.show()

    # Common interface
    def __del__(self):
        d2p.free_all_memory()
        bs.free_all_memory()

    # Interface to spaic2_d2pot
    @property
    def density_params(self):
        """Return c coefficients describing the charge distribtion.
        """
        return d2p.pararho

    @property
    def woodsaxon_params(self):
        """Return B coefficients defining the bumped Wood-Saxon potential.
        """
        return d2p.parapot

    @woodsaxon_params.setter
    def woodsaxon_params(self, params):
        """Set B coefficients defining the bumped Wood-Saxon potential.
        """
        params = np.asarray(params, dtype=np.float64)
        assert len(params.shape) == 1, "received {}".format(params)
        assert all(params >= 0.0) and all(params <= 1.0), "received {}".format(params)
        d2p.set_woodsaxon_params(params)
        # Compute new d2p.pararho aka. density_params
        d2p.compute_potential()
        self.update_density_params()

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
    def quantum_numbers(self):
        return bs.quantum_numbers

    @quantum_numbers.setter
    def quantum_numbers(self, qns):
        qns = np.array(qns, dtype=np.int32)
        prs = np.array(self.density_params, dtype=np.float64)
        is_valid = len(qns.shape) == 2 and qns.shape[1] == 2
        assert is_valid, "QNs example: {} (received {})".format(self.default_quantum_numbers, qns)
        bs.set_density_and_qns(prs, qns)

    def update_density_params(self):
        self.spectra_computed = False
        if self.quantum_numbers is not None:
            qns = np.array(self.quantum_numbers)
            prs = np.array(self.density_params)
            bs.set_density_and_qns(prs, qns)

    def compute_spectra(self):
        bs.compute_beta()
        self.spectra_computed = True
