# Running the TASEP optimization on macOS

The Mac runner activates the correct Conda environment, selects a dataset, limits parallel workers to a Mac-friendly number, prevents the Mac from sleeping, and saves a terminal log.

## One-time setup

From the root of your clone, create the environment if you do not already have it:

```bash
conda env create -f environment.yml
conda activate microlive
```

If `conda activate` is unavailable in a new Terminal window, initialize Conda first:

```bash
source /opt/anaconda3/etc/profile.d/conda.sh
```

## Redraw Figure S11 without running simulations

Open `optimization/plotting_optimization.ipynb` in the analysis environment and run its cells from the repository root or the notebook directory. It reads the four published `results_ACF/results_df_*.csv` grids and writes panel PNGs and SVGs into `optimization_parameters/`. The four `results_ACF/optimized_params_*.csv` files contain the published parameter summaries used downstream for Figure 7 and Table S7.

## Run a new optimization

Starting at the repository root, enter `optimization/` and run:

```bash
cd optimization
./runner_mac.sh DATASET WORKERS
```

`DATASET` can be a name or number:

| Number | Name            | Dataset     |
| -----: | --------------- | ----------- |
|      0 | `utag`        | UTag        |
|      1 | `utag_c_free` | UTag C-free |
|      2 | `suntag`      | SunTag      |
|      3 | `alfatag`     | ALFATag     |

Examples:

```bash
./runner_mac.sh utag
./runner_mac.sh suntag 6
./runner_mac.sh 3 8
```

The worker argument is optional. If omitted, the runner detects the number of logical CPUs and uses a maximum of eight workers. Eight is a reasonable starting point for the current 12-core Mac. Use fewer workers if memory pressure becomes high.

To see the built-in help:

```bash
./runner_mac.sh --help
```

## Input, cache, and output locations

Experimental ACF inputs are loaded from:

```text
optimization/source_files_ACF/data_ACF/
```

Compatible raw simulation caches are reused from, and new simulation caches are written to:

```text
optimization/source_files_ACF/simulation_cache/
```

Final tables and figures are written to:

```text
optimization/results_ACF/
```

Terminal logs are written to:

```text
optimization/results_ACF/logs/
```

Git excludes the simulation cache, logs and rendered plots. The four experimental `df_ACF_*.csv` inputs, their count reports, and the three reporter files in `gene_sequences/utag_project/` are enough to start a new fit with the analysis environment. The original tracking movies are not needed for this step.

The principal outputs are:

- `results_df_<dataset>.csv`: every valid grid-search parameter combination and its cost.
- `optimized_params_<dataset>.csv`: robust parameter estimates and percentile ranges.
- `acf_comp_<dataset>_final.png`: experimental and final simulated ACF comparison.

## Runtime expectations

The complete calculation evaluates an 80 by 80 grid, or 6,400 parameter combinations, using 300 simulation repetitions per combination. A run can therefore take a long time when compatible cached simulations are unavailable. Keep the Mac connected to power and ensure sufficient free storage is available.

The runner uses `caffeinate` so macOS does not sleep during the calculation. Closing Terminal or pressing Control-C stops the active run. Simulations already written to the cache can be reused when the same command is run again.

The final optimized-parameter simulation is deliberately run fresh because `optimization.py` uses `force_fresh=True` for that step. A new stochastic run can differ from the published run and overwrites that dataset's result CSVs. Use a separate checkout for a new fit; use the frozen grids to reproduce the published heatmaps.

## Running with a differently named Conda environment

The default environment is `microlive`. To select another compatible environment:

```bash
UTAG_CONDA_ENV=my_environment ./runner_mac.sh utag 6
```

That environment must contain both `microlive` and `tasep_models`, along with the other packages in `environment.yml`.

