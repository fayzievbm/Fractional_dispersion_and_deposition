"""Reproduce the Chapter 3 Caputo transport and deposition experiments.

Run with Python and NumPy: python chapter3_simulations.py
"""
import csv
from dataclasses import dataclass, replace
from math import gamma
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Parameters:
    length: float = 0.4
    porosity: float = 0.3
    velocity: float = 1e-4
    inlet: float = 0.05
    dispersion: float = 1e-5
    loading: float = 0.005
    capacity: float = 0.1
    initial_rate: float = 32.0
    attachment: float = 42.0
    detachment: float = 0.002


def spatial_matrix(alpha, cells, h):
    """Map full-grid concentrations to the interior Caputo derivative."""
    curvature = np.zeros((cells - 1, cells + 1))
    for i in range(cells - 1):
        curvature[i, i:i + 3] = (1, -2, 1)
    if alpha == 2:
        return curvature / h**2
    k = np.arange(cells)
    weights = np.diff(k**(2 - alpha))
    for i in range(cells - 2, -1, -1):
        curvature[i] = weights[:i + 1] @ curvature[i::-1]
    return curvature / (gamma(3 - alpha) * h**alpha)


def solve(alpha, beta, h, dt, duration=4800, params=Parameters(), ordinary_balance=False):
    """Apply backward Euler to transport and explicit forcing to L1 kinetics."""
    if not (1 < alpha <= 2 and 0 < beta <= 1):
        raise ValueError("Require 1 < alpha <= 2 and 0 < beta <= 1")
    cells, steps = round(params.length / h), round(duration / dt)
    if cells < 2 or not np.isclose(cells * h, params.length) or not np.isclose(steps * dt, duration):
        raise ValueError("The mesh must partition the space and time intervals")
    x, t = np.arange(cells + 1) * h, np.arange(steps + 1) * dt
    space = spatial_matrix(alpha, cells, h)
    operator = params.dispersion * space
    for row in range(cells - 1):
        operator[row, row:row + 2] += np.array([1, -1]) * params.velocity / h
    reduced = operator[:, 1:-1].copy()
    reduced[:, -1] += operator[:, -1]  # C_N = C_(N-1)
    storage = params.porosity / dt
    matrix = storage * np.eye(cells - 1) - reduced
    inverse = np.linalg.solve(matrix, np.eye(cells - 1))
    inlet = operator[:, 0] * params.inlet

    c = np.zeros((steps + 1, cells + 1))
    deposit = np.zeros_like(c)
    increments = np.zeros((steps, cells + 1))
    effective = np.zeros_like(c)
    c[:, 0] = params.inlet
    # Coordinates use reference length 1 m and reference time 1 s.
    history_weights = np.diff(np.arange(steps + 1, dtype=float)**(1 - beta)) if beta < 1 else None
    l1_scale = gamma(2 - beta) * dt**beta
    max_transport_residual = 0.0
    max_kinetic_residual = 0.0
    for n in range(steps):
        mobile, retained = c[n], deposit[n]
        forcing = np.where(
            retained <= params.loading,
            params.initial_rate * params.velocity * mobile,
            params.attachment * params.velocity * mobile - params.detachment * retained,
        )
        forcing[retained >= params.capacity] = 0
        history = history_weights[n:0:-1] @ increments[:n] if beta < 1 and n else 0.0
        increments[n] = l1_scale * forcing - history
        deposit[n + 1] = retained + increments[n]
        l1 = (deposit[n + 1] - retained + history) / l1_scale
        max_kinetic_residual = max(max_kinetic_residual, float(np.max(np.abs(l1 - forcing))))
        sink = increments[n] / dt if ordinary_balance else forcing
        effective[n + 1] = effective[n] + dt * sink
        c[n + 1, 1:-1] = inverse @ (storage * mobile[1:-1] + inlet - sink[1:-1])
        c[n + 1, -1] = c[n + 1, -2]
        residual = storage * (c[n + 1, 1:-1] - mobile[1:-1]) + sink[1:-1] - operator @ c[n + 1]
        max_transport_residual = max(max_transport_residual, float(np.max(np.abs(residual))))
    return dict(x=x, t=t, c=c, deposit=deposit, effective=effective,
                transport_residual=max_transport_residual, kinetic_residual=max_kinetic_residual)


def self_check():
    x = np.linspace(0, 0.4, 41)
    for alpha in (1.6, 1.8, 2):
        operator = spatial_matrix(alpha, 40, 0.01)
        np.testing.assert_allclose(operator @ np.ones_like(x), 0, atol=1e-10)
        np.testing.assert_allclose(operator @ x, 0, atol=1e-10)
        np.testing.assert_allclose(operator @ x**2, 2 * x[1:-1]**(2 - alpha) / gamma(3 - alpha), atol=1e-10)
    for beta in (0.8, 0.9):
        # L1 is exact for a linear history.
        grid = np.arange(21, dtype=float)
        weights = np.diff(grid**(1 - beta))
        approx = np.cumsum(weights) / gamma(2 - beta)
        np.testing.assert_allclose(approx, grid[1:]**(1 - beta) / gamma(2 - beta))
    ordinary = solve(1.8, 1, 0.02, 10, duration=100, ordinary_balance=True)
    fractional = solve(1.8, 1, 0.02, 10, duration=100)
    np.testing.assert_allclose(ordinary["c"], fractional["c"], atol=1e-13)
    np.testing.assert_allclose(ordinary["deposit"], fractional["deposit"], atol=1e-13)


def run():
    self_check()
    output = Path("chapter3_results")
    output.mkdir(exist_ok=True)
    records = []
    cases = [(alpha, beta) for alpha in (2, 1.8, 1.6) for beta in (1, 0.9, 0.8)]
    settings = [
        ("coarse", 0.02, 10.0, Parameters(), False, cases),
        ("medium", 0.01, 5.0, Parameters(), False, cases),
        ("fine", 0.005, 2.5, Parameters(), False, cases),
        ("time", 0.01, 2.5, Parameters(), False, cases),
        ("space", 0.005, 5.0, Parameters(), False, cases),
        ("threshold", 0.005, 2.5, replace(Parameters(), loading=0.008), False, [(2, 0.9), (1.6, 0.8)]),
        ("ordinary", 0.005, 2.5, Parameters(), True, [(2, 0.9), (1.6, 0.8)]),
    ]
    for tag, h, dt, params, ordinary, pairs in settings:
        for alpha, beta in pairs:
            result = solve(alpha, beta, h, dt, params=params, ordinary_balance=ordinary)
            c, rho, x = result["c"], result["deposit"], result["x"]
            assert np.isfinite(c).all() and np.isfinite(rho).all()
            assert result["transport_residual"] < 1e-11 and result["kinetic_residual"] < 1e-11
            assert np.min(c) >= -1e-12 and np.min(rho) >= -1e-12
            assert np.max(c) <= params.inlet + 1e-12
            np.savez_compressed(output / f"{tag}_{alpha:g}_{beta:g}.npz", **result)
            records.append(dict(
                run=tag, alpha=alpha, beta=beta, h=h, dt=dt, loading=params.loading,
                outlet_c=c[-1, -1] / params.inlet,
                rho_002=np.interp(0.02, x, rho[-1]), rho_004=np.interp(0.04, x, rho[-1]),
                rho_020=np.interp(0.20, x, rho[-1]), rho_035=np.interp(0.35, x, rho[-1]),
                inventory=np.trapezoid(rho[-1], x), max_rho=float(rho.max()),
                transport_residual=result["transport_residual"], kinetic_residual=result["kinetic_residual"],
            ))
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    run()
