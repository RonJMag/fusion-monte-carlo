"""Plot the CSV output of fusion_mc in the same style as the original Origin histograms."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data") + os.sep
FIG = os.path.join(ROOT, "figures") + os.sep

BAR = dict(color="#F8CBA0", edgecolor="black", linewidth=0.8)
plt.rcParams.update({"font.size": 11, "axes.labelsize": 12})

def bars(ax, csv, label, xlabel, scale=1.0):
    d = pd.read_csv(csv)
    x = d.iloc[:, 0].values * scale; y = d["count"].values
    w = (x[1] - x[0])
    ax.bar(x, y, width=w, label=label, **BAR)
    ax.set_xlabel(xlabel); ax.set_ylabel("# of Events"); ax.legend(loc="upper right", frameon=True, edgecolor="black")
    return x, y, w

fig, axs = plt.subplots(2, 3, figsize=(17, 10))

# 1. D-T neutron spectrum + Brysk Gaussian
ax = axs[0, 0]
x, y, w = bars(ax, DATA + "DT_thermal_10keV_spectrum_hist.csv", "D-T neutrons, T = 10 keV", r"$K_n$ in MeV", 1e-3)
ev = pd.read_csv(DATA + "DT_thermal_10keV_events.csv")
mu, sd = ev.K3_lab_keV.mean() * 1e-3, ev.K3_lab_keV.std() * 1e-3
brysk_fwhm = 0.5586
g = y.max() * np.exp(-0.5 * ((x - mu) / (brysk_fwhm / 2.3548)) ** 2)
ax.plot(x, g, "k--", lw=1.5, label="Gaussian, Brysk FWHM = 177$\\sqrt{T}$ keV")
ax.legend(loc="upper right", frameon=True, edgecolor="black", fontsize=9)
ax.set_xlim(13.0, 15.1); ax.set_title("Stage 4: Doppler-broadened neutron spectrum")

# 2. spectra at three temperatures
ax = axs[0, 1]
for T, c in [(5, "#F8CBA0"), (10, "#9ECAE1"), (20, "#C7E9C0")]:
    d = pd.read_csv(DATA + f"DT_thermal_{T}keV_spectrum_hist.csv")
    ax.step(d.iloc[:, 0] * 1e-3, d["count"], where="mid", color="black", lw=0.8)
    ax.fill_between(d.iloc[:, 0] * 1e-3, d["count"], step="mid", color=c, alpha=0.8, label=f"T = {T} keV")
ax.set_xlabel(r"$K_n$ in MeV"); ax.set_ylabel("# of Events"); ax.legend(frameon=True, edgecolor="black")
ax.set_xlim(13.0, 15.1); ax.set_title(r"Width $\propto \sqrt{T_i}$: neutron spectrometry as a thermometer")

# 3. reactivity scan
ax = axs[0, 2]
s = pd.read_csv(DATA + "DT_thermal_10keV_reactivity_scan.csv")
ax.errorbar(s.T_keV, s.sv_uniform_m3s, yerr=s.err_uniform, fmt="s", color="#D95F02", label="MC uniform (100000 samples)", capsize=3)
ax.errorbar(s.T_keV, s.sv_importance_m3s, yerr=s.err_importance, fmt="o", mfc="none", color="#1B9E77", label="MC importance-sampled", capsize=3)
Tf = np.logspace(0, 2, 200)
ax.plot(s.T_keV, s.sv_boschhale_m3s, "k-", lw=1.2, label="Bosch-Hale (1992) fit")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"$T_i$ in keV"); ax.set_ylabel(r"$\langle\sigma v\rangle$ in m$^3$ s$^{-1}$")
ax.legend(frameon=True, edgecolor="black", fontsize=9); ax.set_title(r"Stage 1: D-T reactivity $\langle\sigma v\rangle(T)$")

# 4. Gamow peak: reacting E_rel vs Maxwellian
ax = axs[1, 0]
x, y, w = bars(ax, DATA + "DT_thermal_10keV_Erel_hist.csv", "reacting pairs (accepted)", r"$E_{rel}$ in keV")
T = 10.0
fM = np.sqrt(x) * np.exp(-x / T); fM *= y.sum() * w / (fM.sum() * w)
ax.plot(x, fM, "k--", lw=1.5, label="Maxwellian, all pairs")
ax.legend(frameon=True, edgecolor="black", fontsize=9); ax.set_title("Stage 2: Gamow peak from hit-or-miss on $\\sigma v$")

# 5. cos(theta) lab: thermal (flat) vs beam (K_n depends on angle)
ax = axs[1, 1]
b = pd.read_csv(DATA + "DT_beam_100keV_events.csv")
h = ax.hist2d(b.cosTheta3_lab, b.K3_lab_keV * 1e-3, bins=[40, 80], cmap="Oranges")
ax.set_xlabel(r"cosine($\theta_{n,\,lab}$) w.r.t. beam axis"); ax.set_ylabel(r"$K_n$ in MeV")
ax.set_title("Stage 5: 100 keV D beam on 10 keV T plasma")
plt.colorbar(h[3], ax=ax, label="# of Events")

# 6. transport escape spectrum
ax = axs[1, 2]
x, y, w = bars(ax, DATA + "DT_thermal_10keV_transport_escape_hist.csv", "neutrons leaking through 50 cm Li (toy)", r"$K_n$ in MeV")
ax.set_title("Stage 5: toy slab transport, leakage spectrum")

fig.suptitle("Fusion Monte Carlo: 100,000 events per run (generalized from the e$^+$e$^-\\to\\mu^+\\mu^-$ generator)", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(FIG + "fusion_mc_results.png", dpi=150)

# standalone histogram in the exact style of the original two plots
fig, ax = plt.subplots(figsize=(9, 7))
x, y, w = bars(ax, DATA + "DT_thermal_10keV_pT_over_pstar_hist.csv", "Transverse Momenta", r"$p_{T,n}\,/\,p^{*}_{n}$   (lab frame, 10 keV D-T plasma)")
ax.set_title("Neutron transverse momentum: same Jacobian peak as the muon $p_T$ plot")
fig.tight_layout(); fig.savefig(FIG + "fusion_mc_pT_histogram.png", dpi=150)
print("ok")
