def poisson_solver(N_a, N_d, V_bi, x=None, delta=2e-5, n_mesh=50):
    """Numerical solution of the Poisson equation using finite difference."""
    import numpy as np
    q = 1.6e-19; eps_0 = 8.854e-12; eps_r = 11.7; eps = eps_r * eps_0
    if x is None:
        W = np.sqrt((2 * eps * V_bi / q) * (1.0/N_a + 1.0/N_d))
        L = 6.0 * (W / 2.0)
        x = np.linspace(-L, L, n_mesh)
    dx = x[1] - x[0]
    phi = np.zeros_like(x)
    phi[-1] = V_bi
    rho = q * (N_d * 0.5 * (1 + np.tanh(x / delta)) -
               N_a * 0.5 * (1 - np.tanh(x / delta)))
    for _ in range(10000):
        phi_new = phi.copy()
        phi_new[1:-1] = 0.5 * (phi[0:-2] + phi[2:] + dx**2 * rho[1:-1] / eps)
        if np.max(np.abs(phi_new - phi)) < 1e-6:
            break
        phi = phi_new
    E = -np.gradient(phi, dx)
    return E, phi