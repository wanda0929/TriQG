"""
Regenerate the publication-quality pulse-shape figure and the 2-D
fidelity heatmap for the FINAL (strong-blockade, sub-pi/4 area)
configuration:
  Omega_p = 2 pi * 50 MHz, Omega_R = 2 pi * 175 MHz, Omega_c = 2 pi * 50 MHz
  Delta = 2 pi * 500 MHz, alpha = 4, K = 0.95
  sigma = 1.519 ns, T_f = 182.46 ns
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import integrate, optimize

# ---------------------------------------------------------------------
# Pulse model
# ---------------------------------------------------------------------
def pulse(t, T_c, T_f, sigma, A):
    inside = (t >= T_c) & (t <= T_c + 2*T_f)
    tc = T_c + T_f
    expo = -((t - tc)**3 / sigma)**2
    return np.where(inside, (A/2)*np.exp(expo), 0.0)

I_inf = 2 * 0.92770 * 2**(-1.0/6.0)

# Final config
A     = 2*np.pi*50
delta = 2*np.pi*500
T_c   = 0.010
sigma_pi4 = ((2*np.pi*delta)/(A**2 * I_inf))**3
K_final   = 0.95
alpha     = 4.0

sigma_final = K_final**3 * sigma_pi4
T_f_final   = (alpha * sigma_final)**(1.0/3.0)

print(f"sigma_pi4   = {sigma_pi4*1e3:.4f} ns")
print(f"sigma_final = {sigma_final*1e3:.4f} ns,  T_f_final = {T_f_final*1e3:.4f} ns")

# Compare to a few alpha values + the K=1.0 case
configs = [
    # (alpha,        K,    label,                          color)
    (2.4107,       1.00, r"$\alpha=2.41,K=1$ (paper)",     "#8c564b"),
    (4.0,          1.00, r"$\alpha=4,K=1$ (old smooth)",   "#ff7f0e"),
    (4.0,          0.95, r"$\alpha=4,K=0.95$ (FINAL)",     "#2ca02c"),
    (4.0,          0.90, r"$\alpha=4,K=0.90$",             "#1f77b4"),
]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)

# Panel (a): pulse envelopes
ax = axes[0]
for alpha_i, K_i, label, color in configs:
    sigma_i = K_i**3 * sigma_pi4
    T_f_i   = (alpha_i * sigma_i)**(1.0/3.0)
    t = np.linspace(0, 2*T_f_i, 4000)
    f = pulse(t, 0.0, T_f_i, sigma_i, A) / (2*np.pi)
    lw = 2.0 if "FINAL" in label else 1.3
    ax.plot(t*1e3, f, color=color, linewidth=lw, label=label)
ax.set_xlabel("Time t [ns] (target-pulse window, T_c-shifted to 0)")
ax.set_ylabel(r"$\Omega_p(t) / (2\pi)$ [MHz]")
ax.set_title("(a) Super-Gaussian probe envelopes")
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.95)
ax.grid(True, alpha=0.3)

# Panel (b): integrand and cumulative area
ax = axes[1]
ax2 = ax.twinx()
for alpha_i, K_i, label, color in configs:
    sigma_i = K_i**3 * sigma_pi4
    T_f_i   = (alpha_i * sigma_i)**(1.0/3.0)
    t = np.linspace(0, 2*T_f_i, 4000)
    f = pulse(t, 0.0, T_f_i, sigma_i, A)
    integrand = f**2 / (2*delta)
    cum = integrate.cumulative_trapezoid(integrand, t, initial=0)
    lw = 2.0 if "FINAL" in label else 1.3
    ax.plot(t*1e3, integrand, color=color, linewidth=lw, alpha=0.9, label=label)
    ax2.plot(t*1e3, cum, color=color, linewidth=lw, linestyle="--", alpha=0.5)

# pi/4 reference
ax2.axhline(np.pi/4, color="black", lw=0.8, ls=":", alpha=0.6)
ax2.text(0.4, np.pi/4 + 0.01, r"$\pi/4$", fontsize=9, color="black")
ax2.axhline(0.95*np.pi/4, color="#2ca02c", lw=0.8, ls=":", alpha=0.7)

ax.set_xlabel("Time t [ns]")
ax.set_ylabel(r"$\Omega_p^2 / (2\Delta)$  [rad/$\mu$s]    (solid)")
ax2.set_ylabel("Cumulative area    (dashed)")
ax.set_title("(b) Area integrand + cumulative")
ax.legend(loc="upper left", fontsize=8.5)
ax.grid(True, alpha=0.3)

fig.suptitle(
    rf"Probe pulse $\Omega_p(t)$ for the OR gate: $\Omega_p / (2\pi)=50$ MHz, "
    rf"$\Delta / (2\pi)=500$ MHz; final = $\alpha=4$, area $= 0.95\,\pi/4$  "
    rf"($T_f=182.5$ ns, $\sigma=1.52$ ns)",
    fontsize=10.5, y=1.06,
)
fig.savefig("pulse_shapes_final.png", dpi=170, bbox_inches="tight")
print("saved -> pulse_shapes_final.png")
plt.close(fig)

# ---------------------------------------------------------------------
# 2D fidelity heatmap from sweep_2d.log + sweep_2d_smallD.log
# Reuse the runs we already have.
# ---------------------------------------------------------------------
# Hand-typed from the logs (more reliable than re-parsing)
deltas_full = [300, 350, 400, 450, 500, 550, 600, 700, 800, 900]
ks_full     = [0.85, 0.88, 0.90, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.02]
F_data = {
    300: {0.85:0.981130, 0.88:0.986856, 0.90:0.989009, 0.92:0.989877, 0.94:0.989378, 0.96:0.987356, 0.98:0.983644, 1.00:0.978108},
    350: {0.85:0.977789, 0.88:0.985054, 0.90:0.988390, 0.92:0.990435, 0.94:0.991011, 0.96:0.989915, 0.98:0.986988, 1.00:0.982147},
    400: {0.85:0.975004, 0.88:0.983488, 0.90:0.987688, 0.92:0.990535, 0.94:0.991811, 0.96:0.991322, 0.98:0.988933, 1.00:0.984608},
    450: {0.85:0.972637, 0.88:0.982077, 0.90:0.986899, 0.92:0.990312, 0.94:0.992090, 0.96:0.992031, 0.98:0.990053, 1.00:0.986151},
    500: {0.85:0.970602, 0.88:0.980772, 0.90:0.986065, 0.92:0.989899, 0.94:0.992042, 0.95:0.992422, 0.96:0.992324, 0.97:0.991742, 0.98:0.990675, 0.99:0.989133, 1.00:0.987130, 1.02:0.981799},
    550: {0.85:0.968810, 0.88:0.979557, 0.90:0.985213, 0.92:0.989371, 0.94:0.991802, 0.96:0.992353, 0.98:0.990979, 1.00:0.987737},
    600: {0.90:0.984359, 0.93:0.990338, 0.95:0.992066, 0.96:0.992213, 0.97:0.991880, 0.99:0.989804, 1.02:0.983394},
    700: {0.90:0.982693, 0.93:0.989233, 0.95:0.991313, 0.96:0.991638, 0.97:0.991487, 0.99:0.989807, 1.02:0.984086},
    800: {0.90:0.981109, 0.93:0.988050, 0.95:0.990389, 0.96:0.990847, 0.97:0.990836, 0.99:0.989458, 1.02:0.984262},
    900: {0.90:0.979610, 0.93:0.986850, 0.95:0.989387, 0.96:0.989949, 0.97:0.990046, 0.99:0.988903, 1.02:0.984115},
}

# Build matrix (NaN where missing)
M = np.full((len(deltas_full), len(ks_full)), np.nan)
for i, d in enumerate(deltas_full):
    for j, k in enumerate(ks_full):
        if k in F_data[d]:
            M[i, j] = F_data[d][k]

fig, ax = plt.subplots(figsize=(8.5, 4.5), constrained_layout=True)
cmap = plt.cm.viridis.copy()
cmap.set_bad("white")
im = ax.imshow(M, origin="lower", aspect="auto",
               extent=[ks_full[0]-0.01, ks_full[-1]+0.01,
                       deltas_full[0]-25, deltas_full[-1]+25],
               cmap=cmap, vmin=0.97, vmax=0.993)

# Contour at 0.99
KS, DS = np.meshgrid(ks_full, deltas_full)
ax.contour(KS, DS, M, levels=[0.99, 0.992], colors=["white", "yellow"],
           linewidths=[1.2, 1.5], linestyles=["--", "-"])

# Mark the optimum
opt_d, opt_k = 500, 0.95
ax.plot(opt_k, opt_d, "r*", ms=18, mec="white", mew=1.5,
        label=f"optimum: $\\Delta={opt_d}$ MHz, k={opt_k}, "
              f"$\\bar F$=0.9924")
ax.set_xlabel(r"$K = $ area / $(\pi/4)$")
ax.set_ylabel(r"$\Delta / (2\pi)$ [MHz]")
ax.set_title(r"2-D scan of $\bar F_\mathrm{OR}$ over $(\Delta, K)$  "
             r"at $\Omega_p=2\pi\cdot 50$ MHz, $\alpha=4$")
ax.legend(loc="lower right", fontsize=10, framealpha=0.95)

cbar = fig.colorbar(im, ax=ax, label=r"$\bar F_\mathrm{OR}$")
cbar.ax.axhline(0.99, color="white", lw=1.5)

fig.savefig("fidelity_heatmap.png", dpi=170, bbox_inches="tight")
print("saved -> fidelity_heatmap.png")
