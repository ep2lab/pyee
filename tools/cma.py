
import numpy as np
import matplotlib.pyplot as plt
from numpy import logical_and as AND


omega_ce = np.linspace(0, 3, 1000)
omega_pe = np.linspace(0, 3**0.5, 1000)

omega_ce, omega_pe2 = np.meshgrid(np.linspace(0, 3**0.5, 100)**2, np.linspace(0, 3, 100))

cma(omega_pe2, omega_ce, omega_ce, omega_pe2**0.5,0)

def cma(Zm, Rm, omega_ce, omega_pe, lamb):

    omega_pi = omega_pe*np.sqrt(lamb)
    omega_ci = omega_ce*lamb

    R = 1 - omega_pe**2/(1 - omega_ce) - omega_pi**2/(1 + omega_ci)
    L = 1 - omega_pe**2/(1 + omega_ce) - omega_pi**2/(1 - omega_ci)
    P = 1 - omega_pe**2 - omega_pi**2

    # S = 1 - omega_pe**2/(1-omega_ce**2)
    # D = - omega_pe**2*omega_ce/(1 - omega_ce**2)

    S = 1/2*(R + L)
    D = 1/2*(R - L)

    log = np.zeros([Zm.shape[0], Zm.shape[1], 8])
    log[:, :, 6] = np.ones([Zm.shape[0], Zm.shape[1]])
    
    # Region 1
    log[:, :, 0] = R > 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S >= 0
    log[:, :, 4] = D <= 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond1 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))
    # cond1 = D < 0

    # Region 2
    log[:, :, 0] = R < 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D < 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond2 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 3
    log[:, :, 0] = R < 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S < 0
    log[:, :, 4] = D < 0
    log[:, :, 5] = P/S < 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond3 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 4
    log[:, :, 0] = R < 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S < 0
    log[:, :, 4] = D < 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond4 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 5
    log[:, :, 0] = R < 0
    log[:, :, 1] = L < 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S <= 0
    log[:, :, 4] = D <= 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond5 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 6
    log[:, :, 0] = R > 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P > 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D > 0
    log[:, :, 5] = P/S > 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond6 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 7
    log[:, :, 0] = R > 0
    log[:, :, 1] = L > 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D > 0
    log[:, :, 5] = P/S < 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond7 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    # Region 8
    log[:, :, 0] = R > 0
    log[:, :, 1] = L < 0
    log[:, :, 2] = P < 0
    log[:, :, 3] = S > 0
    log[:, :, 4] = D > 0
    log[:, :, 5] = P/S < 0
    # log[:, :, 6] = (R*L - P*S) < 0

    cond8 = AND(log[:, :, 0], AND(log[:, :, 1], AND(log[:, :, 2], AND(log[:, :, 3], AND(log[:, :, 4],AND(log[:, :, 5], log[:, :, 6]))))))

    F = np.zeros(Zm.shape)

    for iz in range(0, Zm[:, 0].size):
        for ix in range(0, Zm[0, :].size):
            if cond1[iz, ix] == 1:
                F[iz, ix] = 1
            elif cond2[iz, ix] == 1:
                F[iz, ix] = 2
            elif cond3[iz, ix] == 1:
                F[iz, ix] = 3
            elif cond4[iz, ix] == 1:
                F[iz, ix] = 4
            elif cond5[iz, ix] == 1:
                F[iz, ix] = 5
            elif cond6[iz, ix] == 1:
                F[iz, ix] = 6
            elif cond7[iz, ix] == 1:
                F[iz, ix] = 7
            elif cond8[iz, ix] == 1:
                F[iz, ix] = 8
            else:
                F[iz, ix] = np.nan

    # scp.gaussian_filter(F , 10)

    p = plt.contourf(Zm, Rm, F, levels = 100)
    plt.colorbar(p, ticks=[1, 2, 3, 4, 5, 6, 7, 8])
    for c in p.collections:
        c.set_edgecolor("r")
    plt.title('Region of propagation')
    plt.show()