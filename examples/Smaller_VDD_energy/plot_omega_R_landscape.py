"""
Plot the (Delta, Omega_R / Omega_p) fidelity landscape from
sweep_omega_R.py + sweep_smallDelta_refine.py.

Two panels:
  (a) Coarse 2D heatmap over (Delta in 200..1000 MHz, ratio in 1..5)
      with F_bar > 0.99 contour overlaid.
  (b) Refinement at Delta = 225 MHz: (ratio, K) heatmap showing the
      knife-edge nature of the small-Delta peak.
"""
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------
# Panel (a): coarse 2D data (hard-coded from sweep_omega_R.log)
# -----------------------------------------------------------------
deltas_A = np.array([200, 300, 400, 500, 700, 1000])
ratios_A = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 5.0])
F_A = np.array([
    [0.863264, 0.983189, 0.904988, 0.991534, 0.978468, 0.979763, 0.852017],
    [0.875265, 0.979391, 0.903543, 0.989293, 0.988740, 0.988568, 0.943285],
    [0.879172, 0.976555, 0.904001, 0.985280, 0.991483, 0.991796, 0.971706],
    [0.880568, 0.974561, 0.904013, 0.982321, 0.991719, 0.992422, 0.982560],
    [0.880809, 0.971668, 0.903069, 0.978281, 0.990292, 0.991313, 0.989011],
    [0.879176, 0.968260, 0.900727, 0.974036, 0.987273, 0.988352, 0.989063],
])

# -----------------------------------------------------------------
# Panel (b): refinement at Delta=225 MHz, (ratio, K) heatmap
# -----------------------------------------------------------------
ratios_C = np.array([2.25, 2.40, 2.50, 2.60, 2.75, 3.00])
Ks_C     = np.array([0.85, 0.90, 0.93, 0.95, 0.97, 1.00])
F_C = np.array([
    [0.898242, 0.930221, 0.947698, 0.957789, 0.966159, 0.974685],
    [0.930959, 0.966366, 0.980893, 0.987080, 0.990252, 0.989262],
    [0.952995, 0.981919, 0.990273, 0.991951, 0.990654, 0.983672],
    [0.967878, 0.986584, 0.988763, 0.986948, 0.982930, 0.973650],
    [0.972801, 0.979527, 0.978121, 0.975748, 0.972578, 0.966674],
    [0.964080, 0.975533, 0.980268, 0.982013, 0.982237, 0.979055],
])

# -----------------------------------------------------------------
# Plot
# -----------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))

# --- Panel (a): coarse heatmap ----------------------------------
ax = axes[0]
levels = np.linspace(0.86, 0.995, 28)
cs = ax.contourf(ratios_A, deltas_A, F_A, levels=levels, cmap="viridis", extend="min")
# overlay F=0.99 and 0.992 contours
ax.contour(ratios_A, deltas_A, F_A, levels=[0.99], colors="white", linewidths=1.6)
ax.contour(ratios_A, deltas_A, F_A, levels=[0.992], colors="red", linewidths=1.3)

# mark the two regimes
ax.plot(3.5, 500, marker="o", mfc="red", mec="white", mew=1.5, ms=10,
        label=r"Farouk: $\bar F=0.9924$, $T_g=385$ ns")
ax.plot(2.5, 200, marker="*", mfc="cyan", mec="black", mew=1.0, ms=14,
        label=r"Small-$\Delta$: $\bar F=0.9915$, $T_g=166$ ns")

ax.set_xlabel(r"$\Omega_R / \Omega_p$")
ax.set_ylabel(r"$\Delta / (2\pi)$  [MHz]")
ax.set_title("(a) OR-gate $\\bar F$ landscape  ($\\Omega_p = 50$ MHz, $K=0.95$)")
ax.legend(loc="upper right", fontsize=9)
cbar = fig.colorbar(cs, ax=ax)
cbar.set_label(r"$\bar F$")

# annotate the ratio=2 dip
ax.text(2.0, 850, r"ratio=2 dip", color="white", fontsize=9, ha="center",
        bbox=dict(facecolor="black", alpha=0.4, edgecolor="none", pad=2))

# --- Panel (b): refinement -------------------------------------
ax = axes[1]
levels2 = np.linspace(0.89, 0.995, 22)
cs2 = ax.contourf(Ks_C, ratios_C, F_C, levels=levels2, cmap="viridis", extend="min")
ax.contour(Ks_C, ratios_C, F_C, levels=[0.99], colors="white", linewidths=1.6)
ax.contour(Ks_C, ratios_C, F_C, levels=[0.992], colors="red", linewidths=1.3)
ax.plot(0.95, 2.5, marker="*", mfc="cyan", mec="black", mew=1.0, ms=14,
        label="best: F = 0.9920")

ax.set_xlabel(r"$K$ = area / $(\pi/4)$")
ax.set_ylabel(r"$\Omega_R / \Omega_p$")
ax.set_title("(b) Robustness at $\\Delta = 225$ MHz  (knife-edge ridge)")
ax.legend(loc="upper right", fontsize=9)
cbar2 = fig.colorbar(cs2, ax=ax)
cbar2.set_label(r"$\bar F$")

# annotate the F>0.99 cells count
ax.text(0.86, 2.95, "Only 4 / 36 cells\nwith $\\bar F > 0.99$",
        color="white", fontsize=9,
        bbox=dict(facecolor="black", alpha=0.55, edgecolor="none", pad=3))

plt.tight_layout()
out = "omega_R_landscape.png"
plt.savefig(out, dpi=160)
print(f"saved -> {out}")
