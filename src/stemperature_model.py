from __future__ import annotations

import numpy as np
from .constants import q, kB, h, m_e_eff, m_h_eff, CM3_TO_M3


def Eg_varshni(T: float) -> float:
    """Bandgap (eV)"""
    return 1.17 - (4.73e-4 * T**2) / (T + 636.0)


def intrinsic_concentration(T: float) -> float:
    """Intrinsic concentration (cm^-3)"""
    Eg_J = Eg_varshni(T) * q

    Nc = 2.0 * ((2*np.pi*m_e_eff*kB*T)/(h**2))**1.5
    Nv = 2.0 * ((2*np.pi*m_h_eff*kB*T)/(h**2))**1.5

    ni_m3 = np.sqrt(Nc * Nv) * np.exp(-Eg_J / (2 * kB * T))
    return ni_m3 / CM3_TO_M3


def built_in_potential(Na: float, Nd: float, T: float) -> float:
    ni = intrinsic_concentration(T)
    return (kB * T / q) * np.log((Na * Nd) / (ni**2))