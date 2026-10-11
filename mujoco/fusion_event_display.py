"""
MuJoCo event display for Fusion Monte Carlo (RonJMag/fusion-monte-carlo).

Reads lab-frame four-vectors from <prefix>_events.csv (output of src/Fusion_Monte_Carlo.cpp)
and launches the reaction products from a point source inside a straight cylindrical
vessel with a uniform axial magnetic field B along z:

  * neutrons (particle 3) are neutral: straight lines until they reach the wall, where
    the hit point is recorded and the particle is stopped ("absorbed").
  * alphas / 3He (particle 4) are charged: the magnetic part of the Lorentz force is applied
    as an exact rotation of v_perp by -omega*dt/2 before and after every MuJoCo step
    (Boris-style split), so they gyrate around field lines with no energy drift.

What MuJoCo does here: integrates the trajectories (rigid point bodies, freejoints, no
gravity, no particle-particle contacts) and renders them. What it does NOT do: nuclear
physics, neutron transport, or self-consistent plasma fields. Those stay in the C++ code.

Time is uniformly rescaled (TIME_SCALE) so a 14 MeV neutron moves a few m/s on screen.
Geometry and B are in real units, and the gyro-frequency is rescaled by the same factor,
so the gyroradii are physical: r_L = p / (|q| B).

Validation checks printed at the end:
  1. Fraction of neutrons hitting the end caps vs the analytic isotropic point-source value
     1 - H / sqrt(H^2 + R^2).
  2. Mean measured alpha gyroradius vs p_perp / (|q| B) from the CSV momenta.
  3. Alpha kinetic-energy drift (a magnetic field does no work, so it should stay ~0).
"""
import argparse, os, json
import numpy as np
import pandas as pd
import mujoco

C = 299792458.0              # m/s
E_CHARGE = 1.602176634e-19   # C
MEV_C_SI = 1e6 * E_CHARGE / C  # 1 MeV/c in kg m/s

ap = argparse.ArgumentParser()
ap.add_argument("--events", default="data/DT_thermal_10keV_events.csv")
ap.add_argument("--n", type=int, default=300, help="number of fusion events to display")
ap.add_argument("--R", type=float, default=2.0, help="vessel radius [m]")
ap.add_argument("--H", type=float, default=2.0, help="vessel half-height [m]")
ap.add_argument("--B", type=float, default=5.0, help="axial field [T]")
ap.add_argument("--Z4", type=int, default=2, help="charge of particle 4 (alpha=2, 3He=2)")
ap.add_argument("--time-scale", type=float, default=1e-7, help="sim speed = real speed * this")
ap.add_argument("--duration", type=float, default=0.8, help="sim seconds")
ap.add_argument("--dt", type=float, default=2.5e-4)
ap.add_argument("--video", default="mujoco/out/fusion_event_display.mp4")
ap.add_argument("--fps", type=int, default=30)
ap.add_argument("--slowmo", type=float, default=6.0, help="video seconds per sim second")
ap.add_argument("--width", type=int, default=960)
ap.add_argument("--height", type=int, default=720)
ap.add_argument("--seed", type=int, default=1)
args = ap.parse_args()
os.makedirs(os.path.dirname(args.video), exist_ok=True)

# ---------------------------------------------------------------- events from the C++ generator
df = pd.read_csv(args.events)
ev = df.sample(n=min(args.n, len(df)), random_state=args.seed).reset_index(drop=True)
p3 = ev[["px3", "py3", "pz3"]].to_numpy() / 1e3     # keV/c -> MeV/c
p4 = ev[["px4", "py4", "pz4"]].to_numpy() / 1e3
E3 = ev["E3"].to_numpy() / 1e3                      # total energy, MeV
E4 = ev["E4"].to_numpy() / 1e3
beta3 = p3 / E3[:, None]                            # v/c = p c / E
beta4 = p4 / E4[:, None]
v3 = beta3 * C * args.time_scale                    # sim m/s
v4 = beta4 * C * args.time_scale
N = len(ev)

# gyro-frequency (relativistic) for particle 4, rescaled to sim time
m4_kg = E4 * MEV_C_SI / C                           # gamma*m = E/c^2  (E in MeV -> kg)
omega4 = args.Z4 * E_CHARGE * args.B / m4_kg * args.time_scale   # rad / sim-s
rL_pred = np.linalg.norm(p4[:, :2], axis=1) * MEV_C_SI / (args.Z4 * E_CHARGE * args.B)

# ---------------------------------------------------------------- MJCF model
R, H = args.R, args.H
bodies = []
for i in range(N):
    bodies.append(f'<body name="n{i}"><freejoint/><geom class="neutron"/></body>')
    bodies.append(f'<body name="a{i}"><freejoint/><geom class="alpha"/></body>')
xml = f"""
<mujoco model="fusion_event_display">
  <size memory="{max(64, N // 10)}M"/>
  <option timestep="{args.dt}" gravity="0 0 0" integrator="Euler">
    <flag contact="disable"/>
  </option>
  <visual>
    <global offwidth="{args.width}" offheight="{args.height}"/>
    <headlight ambient="0.35 0.35 0.35" diffuse="0.6 0.6 0.6"/>
  </visual>
  <default>
    <default class="neutron"><geom type="sphere" size="0.025" mass="1" rgba="0.3 0.7 1 1"/></default>
    <default class="alpha"><geom type="sphere" size="0.03" mass="4" rgba="1 0.55 0.1 1"/></default>
  </default>
  <asset>
    <texture name="sky" type="skybox" builtin="gradient" rgb1="0.02 0.02 0.06" rgb2="0 0 0" width="256" height="256"/>
  </asset>
  <worldbody>
    <light pos="0 0 6" dir="0 0 -1"/>
    <geom name="vessel" type="cylinder" size="{R} {H}" rgba="0.7 0.7 0.8 0.08" contype="0" conaffinity="0"/>
    <geom name="core" type="sphere" size="0.06" rgba="1 0.2 0.6 0.9" contype="0" conaffinity="0"/>
    <geom name="axis" type="cylinder" size="0.005 {H}" rgba="0.9 0.9 0.2 0.4" contype="0" conaffinity="0"/>
    {''.join(bodies)}
  </worldbody>
</mujoco>"""
model = mujoco.MjModel.from_xml_string(xml)
data = mujoco.MjData(model)

nid = np.array([model.body(f"n{i}").id for i in range(N)])
aid = np.array([model.body(f"a{i}").id for i in range(N)])
nq = np.array([model.jnt_qposadr[model.body_jntadr[b]] for b in nid])
aq = np.array([model.jnt_qposadr[model.body_jntadr[b]] for b in aid])
nv_ = np.array([model.jnt_dofadr[model.body_jntadr[b]] for b in nid])
av_ = np.array([model.jnt_dofadr[model.body_jntadr[b]] for b in aid])
mass_a = model.body_mass[aid]

# initial state: all products start at the origin (point plasma source); freejoint qvel[0:3] is world linear velocity
for k in range(N):
    data.qvel[nv_[k]:nv_[k] + 3] = v3[k]
    data.qvel[av_[k]:av_[k] + 3] = v4[k]
mujoco.mj_forward(model, data)

def alpha_vel():
    return np.stack([data.qvel[av_[k]:av_[k] + 3] for k in range(N)])

def pos(qadr):
    return np.stack([data.qpos[q:q + 3] for q in qadr])

def rotate_alphas(theta):
    c, sn = np.cos(theta), np.sin(theta)
    for k in range(N):
        vx, vy = data.qvel[av_[k]], data.qvel[av_[k] + 1]
        data.qvel[av_[k]] = c[k] * vx - sn[k] * vy
        data.qvel[av_[k] + 1] = sn[k] * vx + c[k] * vy

K0 = 0.5 * np.sum(alpha_vel() ** 2, axis=1)
absorbed = np.zeros(N, bool)
hits = np.full((N, 3), np.nan)
a_track_xy = []                                    # for gyroradius measurement
bhat = np.array([0.0, 0.0, 1.0])

renderer = mujoco.Renderer(model, args.height, args.width)
cam = mujoco.MjvCamera()
cam.lookat[:] = [0, 0, 0]; cam.distance = 7.5; cam.elevation = -25; cam.azimuth = 45
frames, next_frame_t = [], 0.0
steps = int(round(args.duration / args.dt))

for s in range(steps):
    # magnetic push for positive charge in B z_hat: dv/dt = omega (v x z_hat) -> rotate v_perp by -omega dt
    rotate_alphas(-0.5 * omega4 * args.dt)
    mujoco.mj_step(model, data)
    rotate_alphas(-0.5 * omega4 * args.dt)

    # neutron wall hits -> record and stop
    pn = pos(nq)
    out = (~absorbed) & ((np.hypot(pn[:, 0], pn[:, 1]) >= R) | (np.abs(pn[:, 2]) >= H))
    for k in np.where(out)[0]:
        hits[k] = pn[k]; absorbed[k] = True
        data.qvel[nv_[k]:nv_[k] + 6] = 0.0
    if s % 4 == 0:
        a_track_xy.append(pos(aq)[:, :2].copy())

    if data.time >= next_frame_t:
        cam.azimuth = 45 + 40 * data.time / args.duration
        renderer.update_scene(data, camera=cam)
        frames.append(renderer.render().copy())
        next_frame_t += 1.0 / (args.fps * args.slowmo)

# ---------------------------------------------------------------- checks
hit_pts = hits[absorbed]
cap = np.abs(hit_pts[:, 2]) >= H - 1e-6
frac_cap = cap.mean() if len(hit_pts) else float("nan")
frac_cap_pred = 1 - H / np.hypot(H, R)
se = np.sqrt(frac_cap_pred * (1 - frac_cap_pred) / max(len(hit_pts), 1))

track = np.stack(a_track_xy)                        # (T, N, 2): circle centre ~ mean, radius ~ mean distance
def fit_circle(xy):          # algebraic least-squares circle fit
    A = np.c_[2 * xy, np.ones(len(xy))]
    b = (xy ** 2).sum(1)
    cx, cy, c0 = np.linalg.lstsq(A, b, rcond=None)[0]
    return np.sqrt(c0 + cx ** 2 + cy ** 2)
rL_meas = np.array([fit_circle(track[:, k, :]) for k in range(N)])
Kf = 0.5 * np.sum(alpha_vel() ** 2, axis=1)
drift = np.abs(Kf / K0 - 1)

summary = {
    "events_file": args.events, "events_displayed": int(N),
    "R_m": R, "H_m": H, "B_T": args.B, "sim_duration_s": args.duration,
    "neutrons_absorbed": int(absorbed.sum()),
    "neutron_endcap_fraction_sim": round(float(frac_cap), 4),
    "neutron_endcap_fraction_isotropic_pred": round(float(frac_cap_pred), 4),
    "endcap_fraction_1sigma": round(float(se), 4),
    "alpha_gyroradius_mean_sim_m": round(float(rL_meas.mean()), 5),
    "alpha_gyroradius_mean_pred_m": round(float(rL_pred.mean()), 5),
    "alpha_gyroradius_ratio": round(float(rL_meas.mean() / rL_pred.mean()), 4),
    "alpha_energy_drift_max": float(drift.max()),
    "frames": len(frames),
    "alpha_gyroradius_max_rel_err": float(np.max(np.abs(rL_meas / rL_pred - 1))),
}
print(json.dumps(summary, indent=2))
with open(os.path.join(os.path.dirname(args.video), "summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
np.savetxt(os.path.join(os.path.dirname(args.video), "neutron_wall_hits.csv"), hit_pts,
           delimiter=",", header="x_m,y_m,z_m", comments="")

import imageio.v2 as imageio
imageio.mimsave(args.video, frames, fps=args.fps, codec="libx264", quality=8)
imageio.imwrite(args.video.replace(".mp4", "_final.png"), frames[-1])
print("wrote", args.video)
