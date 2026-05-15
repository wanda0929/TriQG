"""
Re-verify Rydberg lifetimes for Rb 54 D_{3/2} and Cs 62 D_{5/2} via ARC.

Compares against the n^3-scaled estimates currently used in
`level_parameters.py`:
    tau_R (Rb 54 D_3/2) ~ 74 us   (n^3-scaled from paper's n=66 value)
    tau_r (Cs 62 D_5/2) ~ 77 us   (n^3-scaled from paper's n=76 value)

ARC sources:
    [Sibalic et al., Comp. Phys. Commun. 220, 319 (2017),
     arXiv:1612.05529]
    [Beterov et al., PRA 79, 052504 (2009), arXiv:0902.4995]

We report:
  - 0 K lifetime  (radiative only, spontaneous emission)
  - 300 K lifetime (radiative + BBR-induced)
"""
import arc

print("=" * 72)
print("ARC lifetime re-verification (T = 300 K, BBR included)")
print("=" * 72)
print("ARC version:", arc.__version__)
print()

# ---------- Rb 54 D_{3/2} ----------
rb = arc.Rubidium()
n, l, j = 54, 2, 1.5
tau_Rb_0K   = rb.getStateLifetime(n, l, j, temperature=0,   includeLevelsUpTo=n + 30)
tau_Rb_300K = rb.getStateLifetime(n, l, j, temperature=300, includeLevelsUpTo=n + 30)
print(f"Rb {n} D_{{{int(2*j)}/2}}:")
print(f"    tau (0   K, radiative only) = {tau_Rb_0K * 1e6:7.2f} us")
print(f"    tau (300 K, rad + BBR)      = {tau_Rb_300K * 1e6:7.2f} us")
print(f"    n^3-scaled estimate used     =   74    us")
print(f"    ARC / estimate ratio         = {tau_Rb_300K * 1e6 / 74:.3f}x")
print()

# ---------- Cs 62 D_{5/2} ----------
cs = arc.Caesium()
n, l, j = 62, 2, 2.5
tau_Cs_0K   = cs.getStateLifetime(n, l, j, temperature=0,   includeLevelsUpTo=n + 30)
tau_Cs_300K = cs.getStateLifetime(n, l, j, temperature=300, includeLevelsUpTo=n + 30)
print(f"Cs {n} D_{{{int(2*j)}/2}}:")
print(f"    tau (0   K, radiative only) = {tau_Cs_0K * 1e6:7.2f} us")
print(f"    tau (300 K, rad + BBR)      = {tau_Cs_300K * 1e6:7.2f} us")
print(f"    n^3-scaled estimate used     =   77    us")
print(f"    ARC / estimate ratio         = {tau_Cs_300K * 1e6 / 77:.3f}x")
print()

# ---------- Rb 7 P_{3/2} (intermediate, paper value tau_P = 0.131 us) ----------
n, l, j = 7, 1, 1.5
tau_Rb_P_0K   = rb.getStateLifetime(n, l, j, temperature=0)
tau_Rb_P_300K = rb.getStateLifetime(n, l, j, temperature=300, includeLevelsUpTo=40)
print(f"Rb {n} P_{{{int(2*j)}/2}} (intermediate, low n):")
print(f"    tau (0   K) = {tau_Rb_P_0K * 1e6:7.4f} us")
print(f"    tau (300 K) = {tau_Rb_P_300K * 1e6:7.4f} us")
print(f"    paper value =  0.1310 us")
print()

# ---------- Cross-check the paper's n=66, 76 anchor values ----------
print("-" * 72)
print("Anchor cross-check (Option A levels, paper's n^3 reference):")
print("-" * 72)

# Rb 66 D_{5/2}
n, l, j = 66, 2, 2.5
tau1 = rb.getStateLifetime(n, l, j, temperature=300, includeLevelsUpTo=n + 30)
print(f"Rb 66 D_5/2 @ 300K: {tau1*1e6:6.2f} us   (paper says 134.87 us)")

# Cs 76 D_{3/2}
n, l, j = 76, 2, 1.5
tau2 = cs.getStateLifetime(n, l, j, temperature=300, includeLevelsUpTo=n + 30)
print(f"Cs 76 D_3/2 @ 300K: {tau2*1e6:6.2f} us   (paper says 142.73 us)")
print()
print("=" * 72)
print("Verdict (see report):")
print("=" * 72)
