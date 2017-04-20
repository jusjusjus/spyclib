#!/usr/bin/python3

import spyclib
import numpy as np
import matplotlib.pyplot as plt


# Test spaic2
solver = spyclib.Spaic2Solver()
solver.woodsaxon_params = np.random.rand(30)
solver.plot()

# Test spaic
# solver = spyclib.SpaicSolver()
# solver.plot()
#
# while True:
#     solver.generate_random_potential()
#     solver.plot(show=False)
#     plt.savefig("test.jpg")
#     plt.show()
