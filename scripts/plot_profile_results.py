"""Figure for Stage 6 (tokamak profile source) and the MuJoCo wall-load check."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
pre = "data/DT_profile"
prof = pd.read_csv(f"{pre}_profile.csv"); rh = pd.read_csv(f"{pre}_rho_hist.csv")
sp = pd.read_csv(f"{pre}_spectrum_hist.csv"); summ = pd.read_csv(f"{pre}_profile_summary.csv").set_index("quantity")["value"]
wl = pd.read_csv("mujoco/out/tokamak_stats3000/wall_load_comparison.csv")
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
a = ax[0, 0]; a2 = a.twinx()
a.plot(prof.rho, prof.T_keV, color="#c0392b", lw=2, label="T(ρ) [keV]")
a2.plot(prof.rho, prof.n_m3 / 1e20, color="#2c3e50", lw=2, ls="--", label="n(ρ) [10²⁰ m⁻³]")
a.set_xlabel("normalized minor radius ρ"); a.set_ylabel("T [keV]", color="#c0392b"); a2.set_ylabel("n [10²⁰ m⁻³]", color="#2c3e50")
a.set_title(f"Profiles (P_fus = {float(summ['fusion_power_MW']):.0f} MW, T_burn = {float(summ['T_burn_keV']):.1f} keV)")
a = ax[0, 1]
rg = np.linspace(0, 1, 2001); w = np.interp(rg, prof.rho, prof.emissivity_m3s) * rg
w /= np.trapezoid(w, rg); width = rh.rho.iloc[1] - rh.rho.iloc[0]
a.bar(rh.rho, rh["count"] / rh["count"].sum() / width, width=width, color="#95a5a6", label="sampled birth ρ (100k)")
a.plot(rg, w, color="#e67e22", lw=2, label="ε(ρ)·ρ prediction")
a.set_xlabel("ρ"); a.set_ylabel("probability density"); a.set_title("Where fusion happens"); a.legend()
a = ax[1, 0]
mn, sdp = float(summ["spectrum_mean_keV"]), float(summ["spectrum_sd_pred_keV"])
K = sp.K_keV.to_numpy(); bw = K[1] - K[0]
m3, K30, mD, mT = 939565.42, 14028.0, 1875612.94, 2808921.11
Tg = np.interp(rg, prof.rho, prof.T_keV); sig = np.sqrt(2 * m3 * K30 * Tg / (mD + mT))
mix = np.array([np.trapezoid(w * np.exp(-0.5 * ((k - mn) / sig) ** 2) / (sig * np.sqrt(2 * np.pi)), rg) for k in K])
a.step(K / 1e3, sp["count"] / sp["count"].sum() / bw, where="mid", color="#2c3e50", label="MC volume-integrated")
a.plot(K / 1e3, mix, color="#e67e22", lw=2, label="Brysk mixture over ε-weighted T")
sb = np.sqrt(2 * m3 * K30 * float(summ["T_burn_keV"]) / (mD + mT))
a.plot(K / 1e3, np.exp(-0.5 * ((K - mn) / sb) ** 2) / (sb * np.sqrt(2 * np.pi)), color="#7f8c8d", ls=":", lw=2, label="single Gaussian at T_burn")
a.set_xlim((mn - 1500) / 1e3, (mn + 1500) / 1e3); a.set_xlabel("neutron energy [MeV]"); a.set_ylabel("per keV")
a.set_title("Volume-integrated neutron spectrum"); a.legend(fontsize=8)
a = ax[1, 1]; th = np.degrees(wl.theta_center_rad)
a.plot(th, wl.ref_relative_load, "o-", color="#e67e22", label="ray-cast reference (100k)")
err = wl.mujoco_relative_load / np.sqrt(wl.mujoco_hits)
a.errorbar(th, wl.mujoco_relative_load, yerr=err, fmt="s", color="#2c3e50", capsize=3, label="MuJoCo run (3000)")
a.set_xlabel("poloidal angle on wall [deg]  (0 = outboard midplane)"); a.set_ylabel("relative neutron wall load")
a.set_title("First-wall neutron load"); a.legend()
plt.tight_layout(); plt.savefig("figures/tokamak_profile_results.png", dpi=130); print("ok")
