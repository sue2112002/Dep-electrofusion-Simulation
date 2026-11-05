# Dielectrophoresis + Electrofusion Simulation (MATLAB)

My part of a graduation project: a time-stepped MATLAB simulation of a microfluidic device that **traps** cells with dielectrophoresis (DEP), **pairs** them, and **fuses** pairs with a DC pulse. I used it to test design variants (electrode spacing, where the field is strengthened, cell sizes) before building anything.

> Project name: BioNexus an LOC Device for Rapid Pathogen Detection Integrating AI and Electrochemical Detection  ·  Supervisor: Dr.Sameh Sherif

## What it models
1. **DEP trapping.** Force on each cell `F = 2π·ε_m·r³·Re[CM]·∇(E²)` with a Stokes-drag, time-stepped motion (dt = 10 µs, 0.1 s total). The field along the channel is a cosine profile for interdigitated electrodes, optionally strengthened at chosen electrode pairs.
2. **Pairing.** Trapped cells are paired, giving priority to platelet–red-blood-cell (PLT–RBC) pairs.
3. **Alignment + fusion.** AC alignment phase, then DC pulses decide fused / failed.
4. **Metrics.** Trapping, pairing and fusion efficiency, and the share of PLT–RBC vs same-type fusions.

## Layout
```
src/          dep_ef_simulation_latest.m        most developed version (2nd/4th pairs enhanced +50 %, PLT–RBC priority pairing, efficiency metrics)
              electrode_geometry_field_solver.m finite-difference field for triangular / rectangular / semicircular electrodes (exports CSV)
app/          MicrofluidicElectroporationSystem.m  parameter-panel GUI prototype
experiments/  exp01 … exp08                      design variants (table below)
early_versions/ v1 … v3                          first iterations kept for history
data/         put input_cells_dataset.csv here
```

## Design variants
| File | What changes |
|------|--------------|
| exp01 | Baseline, uniform field, no enhancement |
| exp02 | Field +25 % at electrode pairs 2 and 4 |
| exp03 | Field +150 % at pairs 2 and 4 |
| exp04 | Field +150 % at pairs 3 and 9 |
| exp05 | Asymmetric: +25 % at pair 3, +75 % at pair 9 |
| exp06 | Electrode spacing 40 µm (400 µm channel), +25 % at pairs 2 and 4 |
| exp07 | Spacing 40 µm, lower base field, +100 % at pairs 8 and 9 so only those trap |
| exp08 | Larger cells (synthetic dataset: PLT 1.5–3 µm, RBC 3–7.5 µm) |

## Running
1. MATLAB (the scripts use `readtable` and `xline`, so a recent release).
2. Put `input_cells_dataset.csv` in `data/` (columns: `ID, Type, x, y, radius, sigma_cyt, epsilon_cyt`; `Type` is `PLT` or `RBC`) and run `addpath('data')`. The cell dataset is synthetic (generated for this project), not measured data.
3. Run a script from the `experiments/` or `src/` folder.

## Assumptions and limitations
- Fusion success is an **input probability** (`fusionSuccessRate`, 0.8 in most scripts, 0.95 in the latest), so fusion efficiency reflects that assumption; the model's real outputs are trapping and pairing behaviour.
- The field is an analytic 1-D approximation; Re[CM] uses conductivity only (no frequency dependence).
- The geometry solver builds permittivity and conductivity grids but solves Laplace's equation for a uniform medium.
- Cells that reach a high-field position are snapped to the nearest electrode position.
- In the GUI, the viability tab uses placeholder values and the field tab is illustrative.
