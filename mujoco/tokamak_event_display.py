"""
MuJoCo tokamak event display for Fusion Monte Carlo, Stage 6 (--profile).

Input: <prefix>_events.csv from `fusion_mc --profile`, which carries a birth position
(x, y, z) sampled from the volume emissivity n_D n_T <sigma v>(T(rho)) plus the lab
four-vectors of both products.

Geometry: elongated torus, plasma boundary R = R0 + a cos(t), Z = kappa a sin(t); first wall is
the same shape scaled by --wall (default 1.15 a). Both are drawn as transparent meshes.

Field: B_phi = B0 R0 / R (toroidal), optional poloidal field B_pol = B_phi * r / (q R0) around
the magnetic axis (circular large-aspect-ratio approximation, --q 0 switches it off).

Integration: MuJoCo freejoint point bodies, no gravity, contacts disabled. The magnetic force
is applied as an exact rotation of v about the local b_hat by -omega dt / 2 before and after
each mj_step (Boris-style split), so kinetic energy is conserved to round-off.
Neutrons stop at the first wall; their hit points give the poloidal neutron wall load.

Check (with --q 0 --alphas-only): grad-B + curvature drift in a pure 1/R toroidal field is
vertical with v_d = (v_perp^2 / 2 + v_par^2) / (omega_c R); measured vs predicted.
"""
import argparse, os, json
import numpy as np
import pandas as pd
import mujoco

C = 299792458.0
E_CHARGE = 1.602176634e-19
MEV_KG = 1e6 * E_CHARGE / C**2          # 1 MeV/c^2 in kg

ap = argparse.ArgumentParser()
ap.add_argument("--events", default="data/DT_profile_events.csv")
ap.add_argument("--n", type=int, default=300)
ap.add_argument("--R0", type=float, default=6.2)
ap.add_argument("--a", type=float, default=2.0)
ap.add_argument("--kappa", type=float, default=1.85)
ap.add_argument("--wall", type=float, default=1.15, help="first-wall minor radius / a")
ap.add_argument("--B0", type=float, default=5.3)
ap.add_argument("--q", type=float, default=3.0, help="safety factor for B_pol; 0 = pure toroidal")
ap.add_argument("--Z4", type=int, default=2)
ap.add_argument("--time-scale", type=float, default=1e-7)
ap.add_argument("--duration", type=float, default=1.3)
ap.add_argument("--dt", type=float, default=2.5e-4)
ap.add_argument("--alphas-only", action="store_true", help="skip neutrons (drift check)")
ap.add_argument("--no-video", action="store_true")
ap.add_argument("--outdir", default="mujoco/out/tokamak")
ap.add_argument("--fps", type=int, default=30)
ap.add_argument("--slowmo", type=float, default=5.0)
ap.add_argument("--width", type=int, default=1100)
ap.add_argument("--height", type=int, default=720)
ap.add_argument("--seed", type=int, default=1)
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)
R0, A, KAP, AW = args.R0, args.a, args.kappa, args.wall * args.a

# ---------------------------------------------------------------- events
df = pd.read_csv(args.events)
ev = df.sample(n=min(args.n, len(df)), random_state=args.seed).reset_index(drop=True)
N = len(ev)
x0 = ev[["x_m", "y_m", "z_m"]].to_numpy()
p3 = ev[["px3", "py3", "pz3"]].to_numpy() / 1e3; E3 = ev["E3"].to_numpy() / 1e3   # MeV
p4 = ev[["px4", "py4", "pz4"]].to_numpy() / 1e3; E4 = ev["E4"].to_numpy() / 1e3
v3 = p3 / E3[:, None] * C * args.time_scale
v4 = p4 / E4[:, None] * C * args.time_scale
gm4 = E4 * MEV_KG                                       # gamma m, kg
qm4 = args.Z4 * E_CHARGE / gm4                          # q / (gamma m), C/kg

def bfield(pos):
    """B vector [T] at positions (M,3)."""
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    R = np.hypot(x, y)
    Bphi = args.B0 * R0 / R
    ephi = np.stack([-y / R, x / R, np.zeros_like(R)], 1)
    B = Bphi[:, None] * ephi
    if args.q > 0:
        dR, r = R - R0, np.hypot(R - R0, z)
        r = np.maximum(r, 1e-9)
        Bpol = Bphi * r / (args.q * R0)
        # e_theta (counter-clockwise in the (R, Z) plane): (-z/r) e_R + (dR/r) e_Z
        eR = np.stack([x / R, y / R, np.zeros_like(R)], 1)
        etheta = (-z / r)[:, None] * eR + (dR / r)[:, None] * np.array([0, 0, 1.0])
        B = B + Bpol[:, None] * etheta
    return B

# ---------------------------------------------------------------- meshes (elliptical torus)
def torus_mesh(amin, nt=48, np_=96):
    t = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    ph = np.linspace(0, 2 * np.pi, np_, endpoint=False)
    T, P = np.meshgrid(t, ph, indexing="ij")
    R = R0 + amin * np.cos(T); Z = KAP * amin * np.sin(T)
    V = np.stack([R * np.cos(P), R * np.sin(P), Z], -1).reshape(-1, 3)
    F = []
    for i in range(nt):
        for j in range(np_):
            a_ = i * np_ + j; b = ((i + 1) % nt) * np_ + j
            c = ((i + 1) % nt) * np_ + (j + 1) % np_; d = i * np_ + (j + 1) % np_
            F += [(a_, b, c), (a_, c, d)]
    return " ".join(f"{v:.4f}" for v in V.ravel()), " ".join(str(k) for f in F for k in f)

wv, wf = torus_mesh(AW)
pv, pf = torus_mesh(A)
bodies = []
for i in range(N):
    if not args.alphas_only:
        bodies.append(f'<body name="n{i}" pos="{x0[i,0]} {x0[i,1]} {x0[i,2]}"><freejoint/><geom class="neutron"/></body>')
    bodies.append(f'<body name="a{i}" pos="{x0[i,0]} {x0[i,1]} {x0[i,2]}"><freejoint/><geom class="alpha"/></body>')
xml = f"""
<mujoco model="tokamak_event_display">
  <size memory="{max(64, N // 5)}M"/>
  <option timestep="{args.dt}" gravity="0 0 0" integrator="Euler"><flag contact="disable"/></option>
  <visual><global offwidth="{args.width}" offheight="{args.height}"/>
    <headlight ambient="0.35 0.35 0.35" diffuse="0.6 0.6 0.6"/><map znear="0.01" zfar="200"/></visual>
  <default>
    <default class="neutron"><geom type="sphere" size="0.09" mass="1" rgba="0.3 0.7 1 1"/></default>
    <default class="alpha"><geom type="sphere" size="0.10" mass="4" rgba="1 0.55 0.1 1"/></default>
  </default>
  <asset>
    <mesh name="wall" vertex="{wv}" face="{wf}"/>
    <mesh name="plasma" vertex="{pv}" face="{pf}"/>
  </asset>
  <worldbody>
    <light pos="0 0 20" dir="0 0 -1"/>
    <geom type="mesh" mesh="wall" rgba="0.75 0.75 0.85 0.10" contype="0" conaffinity="0"/>
    <geom type="mesh" mesh="plasma" rgba="1 0.3 0.7 0.07" contype="0" conaffinity="0"/>
    <geom type="cylinder" size="0.15 5" rgba="0.6 0.6 0.6 0.5" contype="0" conaffinity="0"/>
    {''.join(bodies)}
  </worldbody>
</mujoco>"""
model = mujoco.MjModel.from_xml_string(xml)
data = mujoco.MjData(model)

def adr(prefix):
    b = np.array([model.body(f"{prefix}{i}").id for i in range(N)])
    j = model.body_jntadr[b]
    return model.jnt_qposadr[j], model.jnt_dofadr[j]

aq, av = adr("a")
aqi = aq[:, None] + np.arange(3); avi = av[:, None] + np.arange(3)
if not args.alphas_only:
    nq, nv = adr("n")
    nqi = nq[:, None] + np.arange(3); nvi = nv[:, None] + np.arange(3)
    data.qvel[nvi] = v3
data.qvel[avi] = v4
mujoco.mj_forward(model, data)

def rotate_alphas(half):
    """Rotate alpha velocities about local b_hat by angle -omega*dt/2 (positive charge)."""
    pos = data.qpos[aqi]
    B = bfield(pos); Bm = np.linalg.norm(B, axis=1); b = B / Bm[:, None]
    th = -qm4 * Bm * args.time_scale * args.dt * half
    v = data.qvel[avi]
    c, s = np.cos(th)[:, None], np.sin(th)[:, None]
    bv = np.sum(b * v, 1)[:, None]
    # Rodrigues rotation; dv/dt = (q/m) v x B = -omega (b x v), so the angle about b is -omega dt
    data.qvel[avi] = v * c + np.cross(b, v) * s + b * bv * (1 - c)

def inside_wall(pos):
    R = np.hypot(pos[:, 0], pos[:, 1])
    return ((R - R0) / AW) ** 2 + (pos[:, 2] / (KAP * AW)) ** 2 < 1.0

K0 = np.sum(data.qvel[avi] ** 2, 1)
n_alive = np.ones(N, bool) if not args.alphas_only else np.zeros(N, bool)
a_alive = np.ones(N, bool)
n_hits = np.full((N, 3), np.nan)
a_track = []
renderer = None if args.no_video else mujoco.Renderer(model, args.height, args.width)
cam = mujoco.MjvCamera(); cam.lookat[:] = [0, 0, -1.0]; cam.distance = 24; cam.elevation = -38; cam.azimuth = 30
frames, next_t = [], 0.0
steps = int(round(args.duration / args.dt))
for s in range(steps):
    rotate_alphas(0.5)
    mujoco.mj_step(model, data)
    rotate_alphas(0.5)
    if not args.alphas_only:
        pn = data.qpos[nqi]
        out = n_alive & ~inside_wall(pn)
        if out.any():
            n_hits[out] = pn[out]; n_alive[out] = False
            data.qvel[nv[out][:, None] + np.arange(6)] = 0.0
    pa = data.qpos[aqi]
    lost = a_alive & ~inside_wall(pa)
    if lost.any():
        a_alive[lost] = False
        data.qvel[av[lost][:, None] + np.arange(6)] = 0.0
    if s % 8 == 0:
        a_track.append((data.time, pa.copy()))
    if renderer is not None and data.time >= next_t:
        cam.azimuth = 30 + 50 * data.time / args.duration
        renderer.update_scene(data, camera=cam)
        frames.append(renderer.render().copy())
        next_t += 1.0 / (args.fps * args.slowmo)

# ---------------------------------------------------------------- analysis
K1 = np.sum(data.qvel[avi] ** 2, 1)
alive = a_alive
summary = {"events_file": args.events, "events": int(N), "B0_T": args.B0, "q": args.q,
           "R0": R0, "a": A, "kappa": KAP, "wall_minor_radius_m": AW,
           "alpha_energy_drift_max": float(np.max(np.abs(K1[alive] / K0[alive] - 1))) if alive.any() else None,
           "alphas_lost_to_wall": int((~a_alive).sum())}

if not args.alphas_only:
    h = n_hits[~np.isnan(n_hits[:, 0])]
    R = np.hypot(h[:, 0], h[:, 1])
    thw = np.arctan2(h[:, 2] / KAP, R - R0)            # poloidal angle on the wall (0 = outboard midplane)
    nb = 12
    edges = np.linspace(-np.pi, np.pi, nb + 1)
    cnt, _ = np.histogram(thw, edges)
    # wall area per bin: dA = R(t) * dl(t) * 2 pi, dl = AW sqrt(sin^2 + kappa^2 cos^2) dt
    tt = np.linspace(-np.pi, np.pi, 20001); Rt = R0 + AW * np.cos(tt)
    dl = AW * np.sqrt(np.sin(tt) ** 2 + KAP ** 2 * np.cos(tt) ** 2)
    area = np.array([np.trapezoid((2 * np.pi * Rt * dl)[(tt >= edges[k]) & (tt <= edges[k + 1])],
                                  tt[(tt >= edges[k]) & (tt <= edges[k + 1])]) for k in range(nb)])
    load = cnt / area
    centers = 0.5 * (edges[1:] + edges[:-1])
    pd.DataFrame({"theta_center_rad": centers, "hits": cnt, "wall_area_m2": area,
                  "relative_load": load / load.mean()}).to_csv(f"{args.outdir}/neutron_wall_load.csv", index=False)
    pd.DataFrame(h, columns=["x_m", "y_m", "z_m"]).to_csv(f"{args.outdir}/neutron_wall_hits.csv", index=False)
    out_mid = np.abs(centers) < np.pi / 6 + 1e-9
    in_mid = np.abs(np.abs(centers) - np.pi) < np.pi / 6 + 1e-9
    summary.update({"neutrons_hit_wall": int(len(h)),
                    "wall_load_outboard_over_inboard": float(load[out_mid].mean() / load[in_mid].mean())})

if args.q == 0 and len(a_track) > 10:
    # vertical drift check in a pure toroidal field
    t = np.array([tt_ for tt_, _ in a_track]); Z = np.stack([p[:, 2] for _, p in a_track])
    ok = a_alive
    slope = np.array([np.polyfit(t, Z[:, k], 1)[0] for k in range(N)])
    pos0 = x0; Bv = bfield(pos0); Bm = np.linalg.norm(Bv, axis=1); b = Bv / Bm[:, None]
    vpar = np.sum(v4 * b, 1); vperp2 = np.sum(v4 ** 2, 1) - vpar ** 2
    om = qm4 * Bm * args.time_scale
    Rb = np.hypot(pos0[:, 0], pos0[:, 1])
    vd_pred = (0.5 * vperp2 + vpar ** 2) / (om * Rb)
    # direction: b x grad(B)/|grad B| with grad B along -e_R  -> z-component gives the sign (positive charge)
    eR = np.stack([pos0[:, 0] / Rb, pos0[:, 1] / Rb, np.zeros(N)], 1)
    vd_pred = vd_pred * np.cross(b, -eR)[:, 2]
    summary.update({"drift_vz_measured_mean": float(slope[ok].mean()),
                    "drift_vz_predicted_mean": float(vd_pred[ok].mean()),
                    "drift_ratio": float(slope[ok].mean() / vd_pred[ok].mean()),
                    "drift_ratio_per_particle_median": float(np.median(slope[ok] / vd_pred[ok]))})

print(json.dumps(summary, indent=2))
json.dump(summary, open(f"{args.outdir}/summary.json", "w"), indent=2)
if frames:
    import imageio.v2 as imageio
    imageio.mimsave(f"{args.outdir}/tokamak_event_display.mp4", frames, fps=args.fps, codec="libx264", quality=8)
    imageio.imwrite(f"{args.outdir}/tokamak_event_display_final.png", frames[-1])
    imageio.imwrite(f"{args.outdir}/tokamak_event_display_start.png", frames[min(8, len(frames) - 1)])
    print("wrote video", len(frames), "frames")
