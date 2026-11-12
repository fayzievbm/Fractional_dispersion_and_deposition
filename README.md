# Space Fractional Dispersion and Time Fractional Multistage Deposition in Porous Media

Simulation programs, numerical results, and publication figures accompanying our manuscript prepared for *Fractal and Fractional*.

We model suspension transport with a one-sided Caputo spatial derivative and a time-fractional multistage deposition law. The retained phase uses one concentration variable with a piecewise kinetic law for initial loading, attachment–detachment, and a capacity threshold. The repository includes the fractional-balance model and comparisons using an ordinary deposit derivative in the suspension balance.

The published baseline is Fayziev et al., “Modeling the process of anomalous solute transport in a porous medium based on a fractional balance equation,” *WSEAS Transactions on Fluid Mechanics* **20** (2025), 162–168, [doi:10.37394/232013.2025.20.16](https://doi.org/10.37394/232013.2025.20.16). The ordinary-kinetic limit recovers its balance structure. The numerical values here use the present discretization and finite-layer boundary closure.

## Included files

| Location | Contents |
| --- | --- |
| `chapter3_simulations.py` | Main solver and the complete 49-run simulation matrix; the historical filename is retained for compatibility. |
| `make_figures.py` | Generates all 15 manuscript figures from the saved simulation data. |
| `chapter3_results/` | 49 compressed NumPy datasets and `summary.csv`. |
| `figures/` | Figures 1–15 as 4224 × 1488 pixel PNGs with 600 dpi metadata. |
| [downloads/figures-600dpi.zip](downloads/figures-600dpi.zip) | The 15 figures in one submission archive. |
| `verify_results.py` | Checks the stored results, summary, figures, and optional file checksums. |
| `DATA_GUIDE.md` | Parameters, array definitions, result columns, and figure mapping. |
| `legacy/` | The original `Drob_1.py`–`Drob_4.py` and `modul.py`, preserved for provenance. |
| `SHA256SUMS.txt` | Checksums of the published programs, documentation, data, and figures. |

## Run on Windows

1. Download this repository using **Code → Download ZIP**, then extract the ZIP. Alternatively, clone it with Git.
2. Install Python 3.12 or newer if needed. Open PowerShell in the extracted repository folder. All commands below run from that folder.
3. Create a virtual environment and install the dependencies:

   ```powershell
   py -3 -m venv .venv
   .\.venv\Scripts\python.exe -m pip install -r requirements.txt
   ```

   If the `py` launcher is unavailable, use `python -m venv .venv` for the first command. The commands use the environment's Python directly, so activation is unnecessary.

4. Check the downloaded data and figures:

   ```powershell
   .\.venv\Scripts\python.exe verify_results.py --checksums
   ```

5. Regenerate the figures from the included data:

   ```powershell
   .\.venv\Scripts\python.exe make_figures.py
   ```

6. To rerun all simulations, execute the following commands in order:

   ```powershell
   .\.venv\Scripts\python.exe chapter3_simulations.py
   .\.venv\Scripts\python.exe make_figures.py
   .\.venv\Scripts\python.exe verify_results.py
   ```

The full simulation run overwrites the 49 datasets and summary in `chapter3_results/`. It may take several minutes, depending on the computer, and prints no progress messages. Wait until the PowerShell prompt returns before starting plotting. Plotting overwrites the PNGs in `figures/` and prints each saved filename. The submission ZIP is the supplied figure snapshot; rerunning plotting does not update that ZIP.

`--checksums` verifies the exact files as published. After rerunning, numerical validation can pass even if compression, floating-point arithmetic, or fonts produce different file hashes. Run the checker without `--checksums` to validate regenerated data. Windows Arial is used when available; fallback fonts may change typography on other systems.

On Linux or macOS, create an environment with `python3 -m venv .venv` and use `.venv/bin/python` in place of `.\.venv\Scripts\python.exe`.

## Results and checks

The run matrix contains nine spatial/temporal order combinations on five meshes, plus four sensitivity cases. All simulations end at 4800 s. The primary manuscript figures use the fine mesh, with spatial step 0.005 m and time step 2.5 s.

For classical dispersion, reducing the deposition order from 1 to 0.8 changes the final normalized outlet concentration from 0.385919 to 0.032637 and the deposit at 0.04 m from 0.099287 to 0.058321. At deposition order 0.8, reducing the spatial order from 1.8 to 1.6 changes the final outlet concentration from 0.071907 to 0.091037. Full values, mesh comparisons, and sensitivity results are in [summary.csv](chapter3_results/summary.csv).

The solver's built-in checks cover constant, linear, and quadratic spatial profiles, a linear temporal history, and agreement of the two balance formulations at ordinary kinetic order. The stored maximum transport and kinetic residuals over all 49 runs are approximately 3.84 × 10⁻¹⁷ and 3.63 × 10⁻¹⁸. These are discrete-equation checks; the simulations have not been calibrated against experimental retention data. Separate spatial and temporal refinement results quantify numerical sensitivity.

Publication checks were run on Windows with Python 3.12.14, NumPy 2.3.5, and Pillow 13.0.0. The main workflow requires NumPy 2 or newer and Pillow 10.1 or newer. The archived original scripts require additional packages listed in `legacy/requirements.txt` and are not the programs used to generate the current manuscript results.
