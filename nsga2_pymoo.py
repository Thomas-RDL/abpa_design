from matplotlib import pyplot as plt
from pymoo.indicators.kktpm import KKTPM

from battery import *
from pymoo.factory import get_performance_indicator
import time


# NSGA2 PYMOO
#  Created on: jan 2022
#      Author: Thomas Richard de Latour



if __name__ == '__main__':
    mask = ["int", "int", "int"]

    sampling = MixedVariableSampling(mask, {
        "real": get_sampling("real_random"),
        "int": get_sampling("int_random")
    })

    crossover = MixedVariableCrossover(mask, {
        "real": get_crossover("real_sbx", prob=0.9, eta=3.0),
        "int": get_crossover("int_sbx", prob=0.9, eta=3.)
    })

    mutation = MixedVariableMutation(mask, {
        "real": get_mutation("real_pm", eta=3.0),
        "int": get_mutation("int_pm", eta=3.0)
    })

    problem = Battery()
    ngen = 2000
    
    t0 = time.time()
    algorithm = NSGA2(
        pop_size=200,
        n_offsprings=200,
        sampling=sampling,
        crossover=crossover,
        mutation=mutation,
        eliminate_duplicates=True
    )
    termination = get_termination("n_gen", ngen)
    res = minimize(problem,
                   algorithm,
                   termination,
                   seed=1,
                   save_history=True,
                   verbose=True)

    t1 = time.time()

    total = t1 - t0
    print("total time : ", str(total))

    X = res.X
    F = res.F
    H = res.history

    import numpy

    # ******************** Solutions  ****

    a = numpy.asarray(res.X)
    numpy.savetxt("X_nsga_batt.csv", a, delimiter=",")
    a = numpy.asarray(res.F)
    numpy.savetxt("F_nsga_batt.csv", a, delimiter=",")

    Scatter().add(res.F, facecolor="none", edgecolor="red").show()
    Scatter().add(res.F, facecolor="none", edgecolor="red").save("nsga2_batt.png")


