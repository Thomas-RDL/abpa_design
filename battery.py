import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import *
from pymoo.operators.mixed_variable_operator import MixedVariableSampling, MixedVariableMutation, MixedVariableCrossover
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter


# BATTERY PROBLEM PYMOO
#  Created on: jan 2022
#      Author: Thomas Richard de Latour


Eavg = 160  # Wh.km - 1    Avearage    energy    consumption
Tcy = 338.15  # °C    working temperature   of    the    cells
DOD = 50  # % depth of discharge
SOC = 80  # state of charge on average
Crate = 1
PtW = 66.7  # Power to weight ratio (e208)
M_base = 1200  # Mass of the vehicle designed

#    //Ageing coefficients# Qca NMC
a0 = 21500
R = 8.314

# Qca LFP
b0 = 165400
b1 = -4148
b2 = 0.01

# Qcy LFP
d0 = 27789

M_max = 700
V_max = 1
M_carter = 0.55
PtT_efficiency = 90
Int_efficiency = 90
Elec_GWP = 0.0855
LCD = 200000

MANUF_lfp = 8.33  # GWP value for the manufacturing of 1 kg of nmc cell
MANUF_nmc = 17.12  # GWP value for the manufacturing of 1 kg of lfp cell

Ubp_min = 220
Ubp_max = 880

Vef_cyl = 0.295  # Efficiency rate of cell level to system level volumetric cylindrical
Vef_p = 0.353  # Efficiency rate of cell level to system level volumetric prismatic and pouch

Mef_cyl = 0.552  # Efficiency rate of cell level to system level volumetric cylindrical
Mef_p = 0.575  # Efficiency rate of cell level to system level volumetric prismatic and pouch

# Donnée des cellules
name = ["n", "Cc", "Uc", "Tc", "Wc", "Hc", "Mc", "Type", "Form", "Cprice", "Rint", "Ucmax", "Ucmin", "Icm",
        "Crate_peak"]
data = pd.read_csv('cells_data.csv', decimal='.', sep=",")
all_data = pd.DataFrame(data, columns=name)
cell_data = all_data[["Cc", "Uc", "Tc", "Wc", "Hc", "Mc", "Type", "Form", "Cprice", "Rint", "Ucmin", "Crate_peak"]]

class Battery(Problem):
    def __init__(self):
        super().__init__(n_var=3, n_obj=3, n_constr=11)
        self.xl = np.array([0, 1, 1])
        self.xu = np.array([16, 250, 100])

    def _evaluate(self, x, out, *args, **kwargs):
        X = x[:, 0]
        Ns = x[:, 1]
        Np = x[:, 2]

        Ebp = Np * cell_data["Uc"][X] * Ns * cell_data["Cc"][X]
        REV = Ebp * 0.9 / Eavg
        Ncy = LCD / (REV * DOD / 100)

        M_batt = Ns * Np * cell_data["Mc"][X] * cell_data["Form"][X] / Mef_cyl + Ns * Np * cell_data["Mc"][X] * (
                1 - cell_data["Form"][X]) / Mef_p
        E_batt = ((0.3 * LCD * (M_batt * (Ebp / 1000))) / (0.9 * (REV * (M_batt + M_base))))
        E_ri = 0.1 * LCD * (Ebp / 1000) / REV
        GWP = ((0.0855 * E_batt) + (0.0855 * E_ri)) + (
                Ns * Np * cell_data["Mc"][X] * (
                    cell_data["Type"][X] * MANUF_lfp + (1 - cell_data["Type"][X]) * MANUF_nmc))
        Cost = cell_data["Cprice"][X] * Ns * Np
        V_batt = (Ns * Np * cell_data["Hc"][X] * cell_data["Wc"][X] * cell_data["Tc"][X]) * (
                1 - cell_data["Form"][X]) / Vef_p + \
                 Ns * Np * (np.pi * cell_data["Hc"][X] * (cell_data["Wc"][X] ** 2) / 4) * cell_data["Form"][X] / Vef_cyl
        Ppeak = PtW * (M_batt + M_base)

        f1 = GWP
        f2 = -REV
        f3 = Cost

        g0 = 300 - REV
        g1 = 50000 - Ppeak
        g2 = Ppeak - 150000
        g3 = M_batt - M_max
        g4 = Ncy - 1500
        g6 = cell_data["Uc"][X] * Ns - Ubp_max
        g7 = Ubp_min - cell_data["Uc"][X] * Ns
        g8 = Cost - 15000
        g5 = Ppeak - cell_data["Crate_peak"][X] * Ebp + (cell_data["Rint"][X] * (Ns / Np) / 1000) * (
                    (Ppeak / (cell_data["Uc"][X] * Ns)) ** 2)
        g10 = GWP - 7500
        g9 = V_batt - 0.525


        out["F"] = np.column_stack([f1, f2, f3])
        out["G"] = np.column_stack([ g0, g1, g2,g3,g4, g5, g6, g7, g8, g9, g10])

