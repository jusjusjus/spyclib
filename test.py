
import spyclib
import matplotlib.pyplot as plt


solver = spyclib.Spaic2Solver()
solver.plot()


# solver = spyclib.SpaicSolver()
# solver.plot()
#
# while True:
#     solver.generate_random_potential()
#     solver.plot(show=False)
#     plt.savefig("test.jpg")
#     plt.show()
