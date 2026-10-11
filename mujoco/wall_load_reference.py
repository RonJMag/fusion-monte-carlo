"""Independent neutron wall-load reference: straight-line ray casting (NumPy, no MuJoCo) from every
birth point in <prefix>_events.csv along its lab neutron momentum to the elliptical-torus first wall.
Used to check the MuJoCo run's poloidal wall-load histogram."""
import argparse, numpy as np, pandas as pd
ap = argparse.ArgumentParser()
ap.add_argument("--events", default="data/DT_profile_events.csv")
ap.add_argument("--R0", type=float, default=6.2); ap.add_argument("--a", type=float, default=2.0)
ap.add_argument("--kappa", type=float, default=1.85); ap.add_argument("--wall", type=float, default=1.15)
ap.add_argument("--mujoco", default="mujoco/out/tokamak_stats3000/neutron_wall_load.csv")
ap.add_argument("--out", default="mujoco/out/tokamak_stats3000/wall_load_comparison.csv")
a = ap.parse_args(); R0, K, AW = a.R0, a.kappa, a.wall * a.a
df = pd.read_csv(a.events)
x = df[["x_m", "y_m", "z_m"]].to_numpy(); d = df[["px3", "py3", "pz3"]].to_numpy()
d /= np.linalg.norm(d, axis=1)[:, None]
f = lambda p: ((np.hypot(p[:, 0], p[:, 1]) - R0) / AW) ** 2 + (p[:, 2] / (K * AW)) ** 2 - 1.0
lo = np.zeros(len(x)); hi = np.zeros(len(x)); alive = np.ones(len(x), bool); s = 0.0
while alive.any():                                  # coarse march, 2 cm steps, to first exit
    s += 0.02
    out = alive & (f(x + s * d) >= 0)
    hi[out] = s; lo[out] = s - 0.02; alive &= ~out
for _ in range(30):                                 # bisection refine
    m = 0.5 * (lo + hi); o = f(x + m[:, None] * d) >= 0
    hi = np.where(o, m, hi); lo = np.where(o, lo, m)
h = x + hi[:, None] * d
th = np.arctan2(h[:, 2] / K, np.hypot(h[:, 0], h[:, 1]) - R0)
mj = pd.read_csv(a.mujoco); nb = len(mj); edges = np.linspace(-np.pi, np.pi, nb + 1)
cnt, _ = np.histogram(th, edges); load = cnt / mj["wall_area_m2"].to_numpy()
p = cnt / cnt.sum(); n = mj["hits"].sum(); exp = p * n
chi2 = float(((mj["hits"] - exp) ** 2 / exp).sum())
centers = mj["theta_center_rad"].to_numpy()
out_m = np.abs(centers) < np.pi / 6 + 1e-9; in_m = np.abs(np.abs(centers) - np.pi) < np.pi / 6 + 1e-9
res = pd.DataFrame({"theta_center_rad": centers, "ref_relative_load": load / load.mean(),
                    "mujoco_relative_load": mj["relative_load"], "mujoco_hits": mj["hits"], "expected_hits": exp})
res.to_csv(a.out, index=False); print(res.round(3).to_string(index=False))
print(f"reference events {len(x)}; chi2 = {chi2:.1f} for {nb - 1} dof; outboard/inboard ref = {load[out_m].mean() / load[in_m].mean():.3f}")
print(f"max path length {hi.max():.2f} m; mean {hi.mean():.2f} m")
