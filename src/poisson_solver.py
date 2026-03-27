from __future__ import annotations

import numpy as np
from .constants import q, eps_si, CM3_TO_M3, DEFAULT_SMOOTHING
from .analytical_model import depletion_edges


def smooth_step(z, s):
    return 0.5 * (1 + np.tanh(z / s))


def smooth_box(x, a, b, s):
    return smooth_step(x - a, s) - smooth_step(x - b, s)


def charge_density(x, Na, Nd, Vbi, fraction=DEFAULT_SMOOTHING):
    xp, xn, W = depletion_edges(Na, Nd, Vbi)
    s = fraction * W

    Na_m = Na * CM3_TO_M3
    Nd_m = Nd * CM3_TO_M3

    rho_p = -q * Na_m * smooth_box(x, -xp, 0, s)
    rho_n = q * Nd_m * smooth_box(x, 0, xn, s)

    return rho_p + rho_n


def poisson_solver(x, Na, Nd, Vbi, fraction=DEFAULT_SMOOTHING):
    rho = charge_density(x, Na, Nd, Vbi, fraction)

    dx = x[1] - x[0]

    # integrate E
    E = np.zeros_like(x)
    for i in range(1, len(x)):
        E[i] = E[i-1] + 0.5 * (rho[i] + rho[i-1]) / eps_si * dx

    # shift
    E = E - 0.5 * (E[0] + E[-1])

    # integrate phi
    phi = np.zeros_like(x)
    for i in range(1, len(x)):
        phi[i] = phi[i-1] - 0.5 * (E[i] + E[i-1]) * dx

    return E, phi - phi[0]