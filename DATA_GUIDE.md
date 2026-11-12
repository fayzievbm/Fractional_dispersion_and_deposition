# Numerical data guide

## Model and numerical conventions

The computational layer has length 0.4 m, porosity 0.3, velocity 10⁻⁴ m/s, inlet fraction 0.05, and reference dispersion coefficient 10⁻⁵ m²/s. The loading threshold is 0.005 and the nominal capacity is 0.1. Initial-loading and attachment coefficients are 32 and 42; the detachment coefficient is 0.002. Reference length and time are 1 m and 1 s. These reference scales remain fixed when fractional orders change.

The inlet concentration is prescribed, and the outlet concentration equals the adjacent interior node (zero concentration gradient). Initial interior concentration and deposit are zero. The stored inlet concentration is 0.05 at all time indices, including index zero, reflecting the imposed numerical boundary condition.

The spatial order `alpha` is 2, 1.8, or 1.6. The deposition order `beta` is 1, 0.9, or 0.8. Transport uses an implicit update with an upstream advection difference and product-integration spatial weights. Deposition uses an L1 Caputo update with explicit kinetic forcing and the full retained history. The retained state is not clipped at nominal capacity; the threshold switches off instantaneous forcing, while past increments remain in the fractional memory.

In the primary balance, the deposit term is the fractional derivative used in the kinetic equation. In the `ordinary` cases, the deposit term in the suspension balance is the ordinary time difference; fractional kinetics are retained.

## Simulation matrix

Each file is named `{run}_{alpha}_{beta}.npz`.

| Run tag | Spatial step h (m) | Time step dt (s) | Order pairs | Change |
| --- | ---: | ---: | --- | --- |
| `coarse` | 0.020 | 10 | All nine | Coarse mesh |
| `medium` | 0.010 | 5 | All nine | Medium mesh |
| `fine` | 0.005 | 2.5 | All nine | Main figure mesh |
| `time` | 0.010 | 2.5 | All nine | Halve medium time step |
| `space` | 0.005 | 5 | All nine | Halve medium spatial step |
| `threshold` | 0.005 | 2.5 | (2, 0.9), (1.6, 0.8) | Loading threshold increased to 0.008 |
| `ordinary` | 0.005 | 2.5 | (2, 0.9), (1.6, 0.8) | Ordinary deposit term in the balance |

## Arrays inside each NPZ file

Use `numpy.load(path, allow_pickle=False)` to open a dataset.

| Key | Shape | Meaning |
| --- | --- | --- |
| `x` | `(space_nodes,)` | Position in meters, including both endpoints. |
| `t` | `(time_steps + 1,)` | Time in seconds, including zero and 4800 s. |
| `c` | `(time_steps + 1, space_nodes)` | Suspended concentration fraction per liquid volume; divide by 0.05 to obtain the normalized concentration plotted in the figures. |
| `deposit` | Same as `c` | Retained deposit fraction per bulk medium volume. |
| `effective` | Same as `c` | Accumulated `dt * sink` used in the numerical suspension balance. It is a diagnostic of the chosen balance term, not an additional retained compartment. |
| `transport_residual` | Scalar | Maximum absolute discrete transport residual over interior nodes and time steps. |
| `kinetic_residual` | Scalar | Maximum absolute discrete kinetic residual over nodes and time steps. |

For the fractional balance, `effective` is the discrete time sum of the kinetic forcing. For the ordinary balance, it is the accumulated ordinary deposit increment. It should not be interpreted as physical deposit in the fractional cases.

Example:

```python
import numpy as np

with np.load("chapter3_results/fine_1.6_0.8.npz", allow_pickle=False) as data:
    final_outlet = data["c"][-1, -1] / 0.05
    final_deposit = data["deposit"][-1]
    deposit_integral = np.trapezoid(final_deposit, data["x"])
print(final_outlet, deposit_integral)
```

## Summary CSV columns

| Columns | Meaning |
| --- | --- |
| `run`, `alpha`, `beta` | Run tag and fractional orders. |
| `h`, `dt`, `loading` | Spatial step, time step, and loading threshold. |
| `outlet_c` | Final outlet concentration divided by 0.05. |
| `rho_002`, `rho_004`, `rho_020`, `rho_035` | Final physical deposit at 0.02, 0.04, 0.20, and 0.35 m, respectively; interpolated when the location is between grid points. |
| `inventory` | Trapezoidal integral of final physical deposit over the layer, in meters (deposited volume per unit area). |
| `max_rho` | Maximum physical deposit over the complete space–time solution. |
| `transport_residual`, `kinetic_residual` | Maximum absolute discrete equation residuals. |

## Figure mapping

All 15 figures contain two panels, are rendered at 4224 × 1488 pixels, and have 600 dpi metadata.

| Figure | File suffix | Contents |
| ---: | --- | --- |
| 1 | `baseline` | Final profiles at ordinary kinetic order, three spatial orders. |
| 2 | `kinetic` | Final profiles at classical dispersion, three kinetic orders. |
| 3 | `kinetic09` | Classical-dispersion profiles at kinetic order 0.9 and 1600, 3200, 4800 s. |
| 4 | `kinetic08` | Same times at kinetic order 0.8. |
| 5 | `near` | Deposit histories at 0.02 and 0.04 m for the kinetic-order sweep. |
| 6 | `joint1809` | Coupled profiles for (1.8, 0.9) at 450, 900, 1350 s. |
| 7 | `joint1609` | Coupled profiles for (1.6, 0.9), same times. |
| 8 | `joint1808` | Coupled profiles for (1.8, 0.8), same times. |
| 9 | `joint1608` | Coupled profiles for (1.6, 0.8), same times. |
| 10 | `joint` | Final profiles for the four cases with both orders fractional. |
| 11 | `jointnear` | Near-inlet deposit histories for those four cases. |
| 12 | `outlet` | Outlet histories for the kinetic-order sweep and coupled cases. |
| 13 | `inventory` | Physical deposit integrals for those groups. |
| 14 | `variant` | Fractional balance, ordinary balance, and higher threshold for (1.6, 0.8). |
| 15 | `refinement` | Coarse, medium, and fine final profiles for (1.6, 0.8). |

The filenames retain the full form `Figure_01_baseline.png`, etc. The single submission archive is `downloads/figures-600dpi.zip`.
