#!/usr/bin/python3

import spyclib
import numpy as np
import matplotlib.pyplot as plt


# Test cluster radius interface
R = np.array([11.0, 9., 13.])
dR = np.array([0.1, 0.9, 3.])


for i in np.arange(3):
    solver = spyclib.Spaic2Solver() # if this line is moved to the line prior the for-loop, this example does not segfault any more
    print(R[i])
    solver.cluster_radius = R[i]