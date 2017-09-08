#!/usr/bin/python3

import spyclib
import numpy as np
import matplotlib.pyplot as plt


# Test spaic2a
solver = spyclib.Spaic2Solver()
solver.woodsaxon_params = np.random.rand(10)
solver.plot()
