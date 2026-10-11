# MuJoCo event display

`fusion_event_display.py` reads the lab-frame four-vectors written by `src/Fusion_Monte_Carlo.cpp`
(`<prefix>_events.csv`) and moves the reaction products through a straight cylindrical vessel
(radius R, half-height H) with a uniform axial field B, using MuJoCo 3.15 for integration and rendering.

- Neutrons: straight lines to the wall; hit points saved to `neutron_wall_hits.csv`.
- Alphas / 3He: magnetic push applied as an exact v_perp rotation (Boris-style split around each
  `mj_step`), so kinetic energy is conserved to round-off and gyroradii are r_L = p_perp / (|q| B).
- Time is uniformly rescaled (`--time-scale`, default 1e-7); lengths and B are physical.

MuJoCo does not model nuclear physics, transport, or self-consistent fields; those stay in the C++ code.

## Run

    MUJOCO_GL=osmesa PYOPENGL_PLATFORM=osmesa python mujoco/fusion_event_display.py \
        --events data/DT_thermal_10keV_events.csv --n 300 --B 5 --R 2 --H 2

Options: `--events --n --R --H --B --Z4 --time-scale --duration --dt --video --fps --slowmo --seed`.
Use `MUJOCO_GL=egl` or `glfw` where a GPU/display is available.

## Checks (D-T, 10 keV, B = 5 T, R = H = 2 m)

| Quantity | Sim | Prediction |
|---|---|---|
| Neutron end-cap fraction (3000 events) | 0.300 | 1 - H/sqrt(H^2+R^2) = 0.293 ± 0.008 |
| Mean alpha gyroradius (300 events) | 0.04235 m | p_perp/(2eB) = 0.04235 m |
| Max alpha energy drift | 7e-13 | 0 |

# Tokamak event display (Stage 6)

`tokamak_event_display.py` reads `data/DT_profile_events.csv` (from `fusion_mc --profile`). Each
product starts at its sampled birth point inside an ITER-shaped torus (R0 = 6.2 m, a = 2.0 m,
κ = 1.85, first wall at 1.15 a, drawn as transparent meshes).

- Field: B_φ = B0 R0 / R with B0 = 5.3 T, plus an optional poloidal field B_θ = B_φ r / (q R0)
  (large-aspect-ratio circular approximation, default q = 3; `--q 0` turns it off).
- Alphas: an exact Rodrigues rotation of v about the local b̂ by −ω_c dt/2 before and after each
  `mj_step`.
- Neutrons: straight lines to the first wall; the hit points give the poloidal wall load.

    MUJOCO_GL=osmesa PYOPENGL_PLATFORM=osmesa python mujoco/tokamak_event_display.py --n 300 --duration 3.2 --slowmo 3 --width 1104
    MUJOCO_GL=osmesa PYOPENGL_PLATFORM=osmesa python mujoco/tokamak_event_display.py --n 3000 --duration 3.2 --no-video --outdir mujoco/out/tokamak_stats3000
    MUJOCO_GL=osmesa PYOPENGL_PLATFORM=osmesa python mujoco/tokamak_event_display.py --q 0 --alphas-only --no-video --duration 6 --outdir mujoco/out/tokamak_driftcheck
    python mujoco/wall_load_reference.py      # independent NumPy ray-cast check of the wall load

| Check | MuJoCo | Prediction |
|---|---|---|
| Alpha vertical drift in a pure 1/R field (q = 0, 300 alphas, 6 s) | 7.267e-3 m/s (sim units) | (v⊥²/2 + v∥²)/(ω_c R): 7.262e-3, ratio 1.0007 |
| Alpha energy drift, max | ~1e-13 | 0 |
| Poloidal neutron wall load, 12 bins (3000 events) | χ² = 16.6 for 11 dof (p ≈ 0.12) | ray-cast from all 100 000 events |
| Outboard / inboard midplane load | 1.06 (3000) | 1.18 (100k reference) |

Limits: there is no plasma current self-consistency, collisions, or slowing-down, and no neutron
scattering in the wall. The poloidal field is a constant-q approximation.
