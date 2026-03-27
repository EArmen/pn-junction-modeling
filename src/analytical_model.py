from __future__ import annotations

import numpy as np
from .constants import q, eps_si, CM3_TO_M3


def depletion_width(Na: float, Nd: float, Vbi: float) -> float:
    Na_m = Na * CM3_TO_M3
    Nd_m = Nd * CM3_TO_M3

    return np.sqrt((2 * eps_si * Vbi / q) * ((Na_m + Nd_m) / (Na_m * Nd_m)))


def depletion_edges(Na: float, Nd: float, Vbi: float):
    W = depletion_width(Na, Nd, Vbi)
    xp = W * Nd / (Na + Nd)
    xn = W * Na / (Na + Nd)
    return xp, xn, W


def abrupt_field(x, Na, Nd, Vbi):
    Na_m = Na * CM3_TO_M3
    Nd_m = Nd * CM3_TO_M3

    xp, xn, _ = depletion_edges(Na, Nd, Vbi)

    E = np.zeros_like(x)

    mask_p = (x >= -xp) & (x <= 0)
    mask_n = (x > 0) & (x <= xn)

    E[mask_p] = -(q * Na_m / eps_si) * (x[mask_p] + xp)
    E[mask_n] = -(q * Nd_m / eps_si) * (xn - x[mask_n])

    return E


def abrupt_potential(x, Na, Nd, Vbi):
    E = abrupt_field(x, Na, Nd, Vbi)

    phi = np.zeros_like(x)
    for i in range(1, len(x)):
        dx = x[i] - x[i-1]
        phi[i] = phi[i-1] - 0.5 * (E[i] + E[i-1]) * dx

    return phi - phi[0]