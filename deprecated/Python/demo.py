from ttvfaster import run_ttvfaster

params = [
    1.0,  # M_star
    # m1 p1 e1*cos(arg peri1) i1 Omega1 e1*sin(arg peri1) TT1
    0.00001027, 66.03300476, -0.00654273, 1.57079637, -1.57079637, -0.05891280, 142.63816833,
    # m2 p2 e2*cos(arg peri2) i2 Omega2 e2*sin(arg peri2) TT2
    0.00041269, 125.85229492, -0.00096537, 1.57079637, -1.57079637, -0.00953292, 236.66268921
]

print(run_ttvfaster(2, params, 0.0, 1600.0, 6))
