import numpy as np

from ttvfaster import run_ttvfaster

params = np.array([
    1.0,  # M_star
    # m1 p1 e1*cos(arg peri1) i1 Omega1 e1*sin(arg peri1) TT1
    0.00001027, 66.03300476, -0.00654273, 1.57079637, -1.57079637, -0.05891280, 142.63816833,
    0.00041269, 125.85229492, -0.00096537, 1.57079637, -1.57079637, -0.00953292, 236.66268921
])

print(run_ttvfaster(2, params, 0.0, 1600.0, 6))
