
import logging
import numpy as np
from .spaic import spaic

logging.basicConfig(level=logging.INFO)


class SpaicSolver:

    logger = logging.getLogger(name='SpaicSolver')

    default_params = [
            0.189322, 0.776949, 0.284849, 0.706585, 0.780458,
            0.836938, 0.038481, 0.665857, 0.518077, 0.977399
    ]

    def __init__(self):
        self.initialized = False
        self.set_params(self.default_params)

    def __del__(self):
        if self.initialized:
            spaic.free_all_memory()

    def set_params(self, params):
        params = np.asarray(params, dtype=np.float)
        self.params = params
        self.num_params = params.size
        spaic.set_params(self.params)
        if not self.initialized:
            spaic.initialize()
            self.initialized = True
        self.compute_potentials()
        self.frequencies = np.linspace(spaic.wmin, spaic.wmax, spaic.nw)

    def check_asymptotic_points(self, pot):
        condition = abs(pot+(spaic.lin+spaic.dl)*(spaic.lin+spaic.dl+1.0)*spaic.pot1) > 1e-3
        jr = np.arange(spaic.nr)[condition[:-1]][-1]
        assert spaic.nr-jr > 10, "too few asymptotic points: {}".format(spaic.nr-jr)

    def compute_potentials(self):
        wmin = spaic.wmin
        lin  = spaic.lin
        dl   = spaic.dl
        kin  = spaic.kin
        pot = spaic.pot0[:] + np.dot(spaic.pots, 2.0 * self.params - 1.0)
        self.check_asymptotic_points(pot)
        self.bpot = pot[:] + lin * (lin+1.0) * spaic.pot1[:]
        self.cpot = pot[:] + (lin+dl) * (lin+dl+1.0) * spaic.pot1[:]
        self.eb, self.bpsi = spaic.bnumerov(lin, kin, self.bpot)
        assert self.eb+spaic.wmin > 0.0, "Frequency too small!  (eb+wmin={:.3g})".format(self.eb+spaic.wmin)
        assert (np.isnan(self.bpsi) == False).all(), "Computation of bound state failed. (eb={:.3g})".format(self.eb)

    def compute_cpsi(self, E, dl=spaic.dl):
        lout = spaic.lin+dl
        assert not lout < 0, "l_out cannot be negative. (dl = {})".format(dl)
        cpsi = spaic.cnumerov(lout, E, self.cpot, spaic.nr)
        assert (np.isnan(cpsi) == False).all(), "Computation of continuum state failed. (E={:.3g}, dl={})".format(E, dl)
        return cpsi

    def compute_spectrum_at_frequency(self, w, dl=spaic.dl):
        cpsi = self.compute_cpsi(E=self.eb+w, dl=dl)
        return np.sum( self.bpsi * spaic.rr * cpsi )**2

    def compute_spectrum(self):
        spectrum = np.array([
            self.compute_spectrum_at_frequency(w, dl=spaic.dl)
            for w in self.frequencies
        ])
        return spectrum

    def plot(self, **kwargs):
        show = kwargs.get('show', True)
        fs = kwargs.get('fontsize', 14)
        import matplotlib.pyplot as plt
        radius = spaic.rr
        plt.figure(figsize=(15, 7))
        left = 0.05
        bottom = 0.1
        width = 0.65
        height = 1.0 - bottom - 0.05
        # ax = plt.subplot(121)
        ax = plt.axes([left, bottom, width, height])
        plt.title("Potential with Wavefunctions")
        plt.plot(radius, self.bpot, 'k--', lw=1.5, label="Potential")
        plt.plot(radius, self.cpot, 'r--', lw=1.5)
        plt.plot(radius, self.bpsi-0.2, 'b-', lw=2., label="Bound state")
        plt.axhline(y=-0.2, color='b', lw=0.3)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xlabel("Radius", fontsize=fs)
        plt.ylabel("Energie", fontsize=fs)
        E_w = spaic.wmin+0.2*(spaic.wmax-spaic.wmin)
        cpsi = self.compute_cpsi(self.eb+E_w)
        plt.plot(radius, self.eb+E_w + cpsi/10, 'r-', lw=2., label="Unbound state")
        plt.xlim(0, 30)
        plt.legend(loc=0)
        left = left + width + 0.05
        # bottom = 0.1
        width = 1.0 - left - 0.05
        # height = 1.0 - bottom - 0.05
        ax = plt.axes([left, bottom, width, height], sharey=ax, frameon=False, xticks=[], yticks=[])
        plt.title("Spectrum")
        # ax = plt.subplot(122, sharey=ax, frameon=False, xticks=[], yticks=[])
        S = self.compute_spectrum()
        plt.plot(S, self.frequencies, 'k-', lw=2.0)
        # plt.xlim(spaic.wmin, spaic.wmax)
        # plt.ylim(10**-4, None)
        plt.ylim(-0.3, 1.0)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        # plt.xlabel("Frequenz", fontsize=fs)
        plt.setp(ax.get_xticklabels(), visible=False)
        # plt.setp(ax.get_yticklabels(), visible=False)
        # plt.xlabel("Spektrum", fontsize=fs)
        self.write_fortran_config()
        # plt.tight_layout()
        if show:
            plt.show()

    def generate_random_potential(self, num_params=None):
        if num_params is not None:
            self.num_params = num_params
        self.eb = spaic.wmin-1.0
        while True:
            new_params = np.random.rand(self.num_params)
            try: # What if it never works?
                self.set_params(new_params)
                break
            except AssertionError as e:
                self.logger.info('{}, retrying..'.format(e))

    def write_fortran_config(self, filename='./fortran-config.dat'):
        with open(filename, 'w+') as f:
            f.write('{} {} {} {}\n'.format(spaic.rclus, spaic.drclus, spaic.vclus, spaic.fluc))
            f.write('{} {}\n'.format(spaic.lin, spaic.kin))
            f.write('{} {} {} {}\n'.format(spaic.dl, spaic.nw, spaic.wmin, spaic.wmax))
            f.write('{}\n'.format(spaic.npara))
            f.write(' '.join(repr(p) for p in self.params))
