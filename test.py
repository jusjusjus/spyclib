
import spyclib

solver = spyclib.SpaicSolver()
solver.plot()

while True:
    solver.generate_random_potential()
    solver.plot()
