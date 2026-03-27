from __future__ import annotations

# Physical constants
q = 1.602176634e-19          # Coulomb
eps0 = 8.8541878128e-12      # F/m
kB = 1.380649e-23            # J/K
h = 6.62607015e-34           # J*s
m0 = 9.10938356e-31          # kg

# Silicon parameters
epsr_si = 11.7
eps_si = epsr_si * eps0

# Effective masses (simplified model)
m_e_eff = 1.08 * m0
m_h_eff = 0.56 * m0

# Unit conversion
CM3_TO_M3 = 1e6

# Default smoothing parameter
DEFAULT_SMOOTHING = 0.05