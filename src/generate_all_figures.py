import os
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Physical constants
# =========================
q = 1.602176634e-19          # C
eps0 = 8.8541878128e-12      # F/m
epsr_si = 11.7
eps_si = epsr_si * eps0      # F/m
kB = 1.380649e-23            # J/K
h = 6.62607015e-34           # J*s
m0 = 9.10938356e-31          # kg

# Effective masses for a simple silicon classroom model
m_e_eff = 1.08 * m0
m_h_eff = 0.56 * m0

CM3_TO_M3 = 1e6

# =========================
# Utility functions
# =========================
def ensure_plots_dir():
    os.makedirs("plots", exist_ok=True)


def Eg_varshni(T):
    """
    Silicon bandgap (eV), Varshni approximation.
    Gives a realistic decrease of Eg with temperature.
    """
    return 1.17 - (4.73e-4 * T**2) / (T + 636.0)


def intrinsic_concentration(T):
    """
    Intrinsic carrier concentration n_i(T) in cm^-3.
    """
    Eg_J = Eg_varshni(T) * q
    Nc = 2.0 * ((2.0 * np.pi * m_e_eff * kB * T) / (h ** 2)) ** 1.5
    Nv = 2.0 * ((2.0 * np.pi * m_h_eff * kB * T) / (h ** 2)) ** 1.5
    ni_m3 = np.sqrt(Nc * Nv) * np.exp(-Eg_J / (2.0 * kB * T))
    return ni_m3 / CM3_TO_M3


def built_in_potential(Na_cm3, Nd_cm3, T):
    """
    Built-in potential V_bi in volts.
    """
    ni = intrinsic_concentration(T)
    return (kB * T / q) * np.log((Na_cm3 * Nd_cm3) / (ni ** 2))


def depletion_width(Na_cm3, Nd_cm3, Vbi):
    """
    Total depletion width W in meters.
    """
    Na = Na_cm3 * CM3_TO_M3
    Nd = Nd_cm3 * CM3_TO_M3
    return np.sqrt((2.0 * eps_si * Vbi / q) * ((Na + Nd) / (Na * Nd)))


def depletion_edges(Na_cm3, Nd_cm3, Vbi):
    """
    Return x_p, x_n, W in meters.
    """
    W = depletion_width(Na_cm3, Nd_cm3, Vbi)
    xp = W * Nd_cm3 / (Na_cm3 + Nd_cm3)
    xn = W * Na_cm3 / (Na_cm3 + Nd_cm3)
    return xp, xn, W


def abrupt_field_profile(x, Na_cm3, Nd_cm3, Vbi):
    """
    Analytical abrupt-junction electric field profile.
    Returns E(x) in V/m.
    """
    Na = Na_cm3 * CM3_TO_M3
    Nd = Nd_cm3 * CM3_TO_M3
    xp, xn, _ = depletion_edges(Na_cm3, Nd_cm3, Vbi)

    E = np.zeros_like(x)

    mask_p = (x >= -xp) & (x <= 0.0)
    mask_n = (x > 0.0) & (x <= xn)

    # Physical field sign convention
    E[mask_p] = -(q * Na / eps_si) * (x[mask_p] + xp)
    E[mask_n] = -(q * Nd / eps_si) * (xn - x[mask_n])

    return E


def abrupt_potential_profile(x, Na_cm3, Nd_cm3, Vbi):
    """
    Analytical abrupt-junction potential profile.
    Returns phi(x) in volts.
    """
    E = abrupt_field_profile(x, Na_cm3, Nd_cm3, Vbi)

    phi = np.zeros_like(x)
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        phi[i] = phi[i - 1] - 0.5 * (E[i] + E[i - 1]) * dx

    phi = phi - phi[0]
    return phi


def smooth_step(z, s):
    return 0.5 * (1.0 + np.tanh(z / s))


def smooth_box(x, a, b, s):
    """
    Smooth approximation of a box function equal to 1 on [a,b].
    """
    return smooth_step(x - a, s) - smooth_step(x - b, s)


def smooth_charge_density(x, Na_cm3, Nd_cm3, Vbi, fraction=0.05):
    """
    Smoothed finite-support charge profile rho(x) in C/m^3.
    The smoothing width is delta = fraction * W.
    """
    xp, xn, W = depletion_edges(Na_cm3, Nd_cm3, Vbi)
    s = fraction * W

    Na = Na_cm3 * CM3_TO_M3
    Nd = Nd_cm3 * CM3_TO_M3

    rho_p = -q * Na * smooth_box(x, -xp, 0.0, s)
    rho_n = +q * Nd * smooth_box(x, 0.0, xn, s)

    return rho_p + rho_n


def numerical_profiles(x, Na_cm3, Nd_cm3, Vbi, fraction=0.05):
    """
    Compute smoothed numerical E(x), phi(x) by integrating rho.
    Returns:
        E_num (V/m)
        phi_num (V)
    """
    rho = smooth_charge_density(x, Na_cm3, Nd_cm3, Vbi, fraction=fraction)

    dx = x[1] - x[0]

    # Integrate Poisson: dE/dx = rho/eps
    E = np.zeros_like(x)
    for i in range(1, len(x)):
        E[i] = E[i - 1] + 0.5 * (rho[i] + rho[i - 1]) / eps_si * dx

    # Shift so that the field is approximately centered around zero outside depletion
    E = E - 0.5 * (E[0] + E[-1])

    # Integrate E = -dphi/dx
    phi = np.zeros_like(x)
    for i in range(1, len(x)):
        phi[i] = phi[i - 1] - 0.5 * (E[i] + E[i - 1]) * dx

    phi = phi - phi[0]

    return E, phi


def eps_max(y_num, y_ref):
    return np.max(np.abs(y_num - y_ref))


def eps_rms(y_num, y_ref):
    return np.sqrt(np.mean((y_num - y_ref) ** 2))


# =========================
# Figure 1a, 1b
# =========================
def make_fig1():
    T_values = np.linspace(200, 400, 300)
    Na = 1e18
    Nd = 1e18

    ni_values = np.array([intrinsic_concentration(T) for T in T_values])
    vbi_values = np.array([built_in_potential(Na, Nd, T) for T in T_values])

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(T_values, vbi_values, linewidth=2)
    ax.set_xlabel("Temperature (K)", fontsize=16)
    ax.set_ylabel(r"Built-in Potential, $V_{bi}$ (V)", fontsize=16)
    ax.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/1a.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(T_values, ni_values, linewidth=2)
    ax.set_yscale("log")
    ax.set_xlabel("Temperature (K)", fontsize=16)
    ax.set_ylabel(r"Intrinsic Carrier Concentration, $n_i$ (cm$^{-3}$)", fontsize=16)
    ax.grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/1b.png", dpi=200)
    plt.close(fig)


# =========================
# Figure 2
# =========================
def make_fig2():
    doping = np.logspace(16, 20, 300)
    T = 300.0

    widths = []
    for N in doping:
        Vbi = built_in_potential(N, N, T)
        W = depletion_width(N, N, Vbi)
        widths.append(W)

    widths = np.array(widths)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(doping, widths, linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel(r"Doping Concentration, $N_a = N_d$ (cm$^{-3}$)", fontsize=16)
    ax.set_ylabel(r"Total Depletion Width, $W$ (m)", fontsize=16)
    ax.grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/2.png", dpi=200)
    plt.close(fig)


# =========================
# Figure 3a, 3b
# =========================
def make_fig3():
    T = 300.0
    cases = [1e16, 1e18, 1e20]
    colors = ["tab:blue", "tab:green", "tab:red"]

    # Use a common x-axis in micrometers
    x_um = np.linspace(-1.4, 1.4, 3000)
    x_m = x_um * 1e-6

    figE, axE = plt.subplots(figsize=(10, 7))
    figP, axP = plt.subplots(figsize=(10, 7))

    for N, c in zip(cases, colors):
        Vbi = built_in_potential(N, N, T)
        E = abrupt_field_profile(x_m, N, N, Vbi) / 1e3  # kV/m
        phi = abrupt_potential_profile(x_m, N, N, Vbi)

        axE.plot(x_um, E, color=c, linewidth=2,
                 label=rf"$N_a = N_d = 10^{{{int(np.log10(N))}}}\ \mathrm{{cm}}^{{-3}}$")

        axP.plot(x_um, phi, color=c, linewidth=2,
                 label=rf"$N_a = N_d = 10^{{{int(np.log10(N))}}}\ \mathrm{{cm}}^{{-3}}$")

    axE.set_xlabel(r"Position, $x$ ($\mu$m)", fontsize=16)
    axE.set_ylabel(r"Electric Field, $E$ (kV/m)", fontsize=16)
    axE.legend(fontsize=14)
    axE.grid(True, alpha=0.35)

    axP.set_xlabel(r"Position, $x$ ($\mu$m)", fontsize=16)
    axP.set_ylabel(r"Potential, $\phi$ (V)", fontsize=16)
    axP.legend(fontsize=14)
    axP.grid(True, alpha=0.35)

    plt.figure(figE.number)
    plt.tight_layout()
    plt.savefig("plots/3a.png", dpi=200)
    plt.close(figE)

    plt.figure(figP.number)
    plt.tight_layout()
    plt.savefig("plots/3b.png", dpi=200)
    plt.close(figP)


# =========================
# Figure 4a, 4b
# =========================
def make_fig4():
    T = 300.0
    Na = 1e18
    Nd = 1e16
    Vbi = built_in_potential(Na, Nd, T)

    xp, xn, W = depletion_edges(Na, Nd, Vbi)
    L = 1.1 * W

    x_m = np.linspace(-L, L, 3000)
    x_um = x_m * 1e6

    E_num, phi_num = numerical_profiles(x_m, Na, Nd, Vbi, fraction=0.05)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(x_um, E_num / 1e3, linewidth=2)
    ax.set_xlabel(r"Position, $x$ ($\mu$m)", fontsize=16)
    ax.set_ylabel(r"Electric Field, $E$ (kV/m)", fontsize=16)
    ax.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/4a.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(x_um, phi_num, linewidth=2)
    ax.set_xlabel(r"Position, $x$ ($\mu$m)", fontsize=16)
    ax.set_ylabel(r"Potential, $\phi$ (V)", fontsize=16)
    ax.grid(True, alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/4b.png", dpi=200)
    plt.close(fig)


# =========================
# Figure 5
# =========================
def make_fig5():
    T = 300.0
    cases = [
        (r"$10^{16}$", 1e16, 1e16),
        (r"$10^{18}$", 1e18, 1e18),
        (r"$10^{20}$", 1e20, 1e20),
        (r"$10^{18},10^{16}$", 1e18, 1e16),
    ]

    emax_values = []
    erms_values = []
    labels = []

    for label, Na, Nd in cases:
        Vbi = built_in_potential(Na, Nd, T)
        xp, xn, W = depletion_edges(Na, Nd, Vbi)
        L = 1.2 * W

        x = np.linspace(-L, L, 3000)
        E_ref = abrupt_field_profile(x, Na, Nd, Vbi)
        E_num, _ = numerical_profiles(x, Na, Nd, Vbi, fraction=0.05)

        labels.append(label)
        emax_values.append(eps_max(E_num, E_ref))
        erms_values.append(eps_rms(E_num, E_ref))

    x_pos = np.arange(len(labels))
    width = 0.36

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(x_pos - width / 2, emax_values, width=width, color="blue", label=r"$\varepsilon_{\max}$")
    ax.bar(x_pos + width / 2, erms_values, width=width, color="red", label=r"$\varepsilon_{\mathrm{RMS}}$")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel(r"$N_a, N_d$ (cm$^{-3}$)", fontsize=16)
    ax.set_ylabel(r"Error (V/m)", fontsize=16)
    ax.legend(fontsize=14)
    ax.grid(True, axis="y", alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/5.png", dpi=200)
    plt.close(fig)


# =========================
# Figure 6
# =========================
def make_fig6():
    T = 300.0
    Na = 1e18
    Nd = 1e18
    Vbi = built_in_potential(Na, Nd, T)
    xp, xn, W = depletion_edges(Na, Nd, Vbi)

    x = np.linspace(-1.2 * W, 1.2 * W, 3000)
    E_ref = abrupt_field_profile(x, Na, Nd, Vbi)

    fractions = np.array([0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70])
    errors = []

    for f in fractions:
        E_num, _ = numerical_profiles(x, Na, Nd, Vbi, fraction=f)
        errors.append(eps_rms(E_num, E_ref))

    errors = np.array(errors)

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.plot(fractions, errors, marker="o", linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel(r"Smoothing fraction, $f$ in $\delta = f \cdot W$", fontsize=16)
    ax.set_ylabel(r"RMS Error vs Reference, $\varepsilon_{\mathrm{RMS}}$ (V/m)", fontsize=16)
    ax.grid(True, which="both", alpha=0.35)
    plt.tight_layout()
    plt.savefig("plots/6.png", dpi=200)
    plt.close(fig)


# =========================
# Main
# =========================
def main():
    ensure_plots_dir()
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    make_fig5()
    make_fig6()
    print("All figures have been generated in the folder: plots/")


if __name__ == "__main__":
    main()