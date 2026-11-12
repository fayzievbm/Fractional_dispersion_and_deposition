"""Validate the supplied numerical data and publication figures.

Run from the repository folder: python verify_results.py [--checksums]
"""
import argparse
import csv
import hashlib
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from PIL import Image

from chapter3_simulations import self_check


ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checksums', action='store_true',
                        help='Also verify the exact published files.')
    args = parser.parse_args()
    self_check()
    with (ROOT / 'chapter3_results/summary.csv').open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    pairs = {(a, b) for a in (2, 1.8, 1.6) for b in (1, 0.9, 0.8)}
    expected = {(tag, a, b) for tag in ('coarse', 'medium', 'fine', 'time', 'space')
                for a, b in pairs}
    expected |= {(tag, a, b) for tag in ('threshold', 'ordinary')
                 for a, b in ((2, 0.9), (1.6, 0.8))}
    found = {(r['run'], float(r['alpha']), float(r['beta'])) for r in rows}
    require(len(rows) == 49 and found == expected, 'Unexpected simulation matrix')
    npz_paths = set((ROOT / 'chapter3_results').glob('*.npz'))
    expected_paths = set()
    worst_transport = worst_kinetic = 0.0
    for row in rows:
        a, b = float(row['alpha']), float(row['beta'])
        path = ROOT / 'chapter3_results' / f'{row["run"]}_{a:g}_{b:g}.npz'
        expected_paths.add(path)
        with np.load(path, allow_pickle=False) as data:
            x, t, c, rho = (data[k] for k in ('x', 't', 'c', 'deposit'))
            h, dt = float(row['h']), float(row['dt'])
            np.testing.assert_allclose(x, np.arange(round(0.4 / h) + 1) * h, atol=1e-13)
            np.testing.assert_allclose(t, np.arange(round(4800 / dt) + 1) * dt, atol=1e-13)
            for key in ('c', 'deposit', 'effective'):
                require(data[key].shape == (len(t), len(x)), f'{path.name}: {key} shape')
                require(np.isfinite(data[key]).all(), f'{path.name}: non-finite {key}')
            require(c.min() >= -1e-12 and c.max() <= 0.05 + 1e-12,
                    f'{path.name}: concentration bounds')
            require(rho.min() >= -1e-12, f'{path.name}: negative deposit')
            np.testing.assert_allclose(c[:, 0], 0.05, atol=1e-13)
            np.testing.assert_allclose(c[:, -1], c[:, -2], atol=1e-13)
            np.testing.assert_allclose(rho[0], 0, atol=1e-13)
            values = dict(outlet_c=c[-1, -1] / 0.05,
                          inventory=np.trapezoid(rho[-1], x), max_rho=rho.max())
            for field, point in (('rho_002', 0.02), ('rho_004', 0.04),
                                 ('rho_020', 0.20), ('rho_035', 0.35)):
                values[field] = np.interp(point, x, rho[-1])
            for field in ('transport_residual', 'kinetic_residual'):
                values[field] = float(data[field])
                require(0 <= values[field] < 1e-11, f'{path.name}: {field}')
            for field, value in values.items():
                np.testing.assert_allclose(value, float(row[field]), rtol=1e-10, atol=1e-13,
                                           err_msg=f'{path.name}: {field}')
            worst_transport = max(worst_transport, values['transport_residual'])
            worst_kinetic = max(worst_kinetic, values['kinetic_residual'])
    require(npz_paths == expected_paths, 'Missing or extra NPZ files')
    figures = sorted((ROOT / 'figures').glob('Figure_*.png'))
    require(len(figures) == 15, 'Expected 15 figures')
    for i, path in enumerate(figures, 1):
        require(path.name.startswith(f'Figure_{i:02d}_'), 'Figure numbering gap')
        with Image.open(path) as im:
            require(im.size == (4224, 1488), f'{path.name}: image dimensions')
            require(min(im.info.get('dpi', (0, 0))) >= 599.9, f'{path.name}: DPI metadata')
            im.verify()
    with ZipFile(ROOT / 'downloads/figures-600dpi.zip') as archive:
        require(archive.testzip() is None, 'Figure archive is corrupt')
        require(set(archive.namelist()) == {p.name for p in figures}, 'Figure archive contents')
    if args.checksums:
        lines = (ROOT / 'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines()
        for line in lines:
            expected_hash, name = line.split('  ', 1)
            path = (ROOT / name).resolve()
            require(path.is_relative_to(ROOT), 'Checksum path outside repository')
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            require(actual == expected_hash, f'Checksum mismatch: {name}')
        print(f'PASS: {len(lines)} file checksums')
    print('PASS: solver self-checks; 49 datasets and summary; 15 figures and submission ZIP')
    print(f'Maximum stored residuals: transport={worst_transport:.3e}, kinetic={worst_kinetic:.3e}')


if __name__ == '__main__':
    main()
