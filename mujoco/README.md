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
