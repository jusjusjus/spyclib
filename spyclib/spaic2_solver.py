
import logging
import numpy as np
from .spaic2 import spaic2_betaspect as bs
from .spaic2 import spaic2_d2pot as d2p

logging.basicConfig(level=logging.INFO)
EPS = 1.0**-6


class Spaic2Solver:

    logger = logging.getLogger(name='Spaic2Solver')

    default_density_params = np.array([
            0.84, 0.39, 0.78, 0.80, 0.91, 0.20, 0.34, 0.77, 0.28,
            0.55, 0.48, 0.63, 0.36, 0.51, 0.95, 0.92, 0.64, 0.72, 0.14, 0.61,
            0.02, 0.24, 0.14, 0.80, 0.16, 0.40, 0.13, 0.11, 1.00, 0.22
    ])

    default_quantum_numbers = np.array([
            [0, 1],
            [1, 1]
    ], dtype=np.int32)

    def __init__(self):
        self.density_params   = self.default_density_params # setter
        self.quantum_numbers  = self.default_quantum_numbers # setter

    def plot(self, **kwargs):
        if not self.spectra_computed: self.compute_spectra()

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
    def potential_params(self):
        return d2p.pararho

    @property
    def density_params(self):
        return d2p.parapot

    @density_params.setter
    def density_params(self, params):
        params = np.asarray(params, dtype=np.float64)
        assert len(params.shape) == 1, "recieved {}".format(params)
        assert all(params >= 0.0) and all(params <= 1.0), "recieved {}".format(params)
        d2p.set_density_params(params)
        # Compute new d2p.pararho aka. potential_params
        d2p.compute_potential()
        self.update_potential_params()

    # Interface to spaic2_betaspect
    @property
    def beta_spectra(self):
        return bs.beta

    @property
    def spectra(self):
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
        prs = np.array(self.potential_params, dtype=np.float64)
        is_valid = len(qns.shape) == 2 and qns.shape[1] == 2
        assert is_valid, "QNs example: {} (received {})".format(self.default_quantum_numbers, qns)
        bs.set_potential_and_qns(prs, qns)

    def update_potential_params(self):
        self.spectra_computed = False
        if self.quantum_numbers is not None:
            qns = np.array(self.quantum_numbers)
            prs = np.array(self.potential_params)
            bs.set_potential_and_qns(prs, qns)

    def compute_spectra(self):
        bs.compute_beta()
        self.spectra_computed = True
