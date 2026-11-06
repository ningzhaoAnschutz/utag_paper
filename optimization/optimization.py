#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Estimate translation initiation and elongation rates from experimental ACFs.

The script evaluates a grid of initiation (ki) and elongation (ke) rates with
TASEP simulations, compares each simulated autocorrelation function (ACF) with
the selected experimental ACF, and saves the resulting parameter landscape.

Usage:
    python optimization.py [dataset_index]

Dataset indices:
    0: UTag
    1: UTag cysteine-free
    2: SunTag
    3: ALFATag
"""

import hashlib
import json
import os
import pickle
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

matplotlib.use("Agg")  # Use a non-interactive backend on clusters and macOS.
import matplotlib.pyplot as plt

# Plot styling
plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.linewidth": 1.5,
        "font.family": "sans-serif",
        "font.sans-serif": "Arial",
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.labelcolor": "black",
        "text.color": "black",
        "xtick.color": "black",
        "ytick.color": "black",
    }
)

# Dataset selection and analysis settings
DATASET_NAMES = {
    0: "utag",
    1: "utag_c_free",
    2: "suntag",
    3: "alfatag",
}

dataset_index = int(sys.argv[1]) if len(sys.argv) > 1 else 0
if dataset_index not in DATASET_NAMES:
    raise ValueError(f"Invalid dataset index: {dataset_index}. Must be 0-3.")
dataset = DATASET_NAMES[dataset_index]

STEP_SIZE_IN_SEC = 5
DOWNSAMPLE = False
DOWNSAMPLING_FACTOR = 2
Y_AXIS_LIMITS = [-0.01, 0.1]
X_AXIS_LIMITS = [-0.2, 25.5]
MAD_THRESHOLD_FACTOR = 4
NUMBER_REPETITIONS = 300
BURN_IN_TIME_SEC = 2000
SIMULATION_DURATION_SEC = 360 * 5
GRID_SIZE = 80
KI_MIN = 0.01
KI_MAX = 0.08
KE_MIN = 2
KE_MAX = 10
TOP_PERCENTILE = 5.0
CONFIDENCE_LEVEL = 0.90
PARALLEL_JOBS = int(os.environ.get("UTAG_N_JOBS", "88"))

if PARALLEL_JOBS == 0 or PARALLEL_JOBS < -1:
    raise ValueError("UTAG_N_JOBS must be -1 or a positive integer.")


def find_project_root():
    """Find the UTag project root directory."""
    current_dir = Path(__file__).resolve().parent
    check_dir = current_dir

    while check_dir != check_dir.parent:
        if check_dir.name in ("utag", "utag_paper") and (check_dir / "src").exists():
            return check_dir
        if (check_dir / "src").exists() and (check_dir / "optimization").exists():
            return check_dir
        check_dir = check_dir.parent

    possible_locations = [
        Path.home() / "Desktop" / "utag_paper",
        Path.home() / "Desktop" / "utag",
        Path.home() / "utag_paper",
        Path.home() / "utag",
        Path.cwd().parent if Path.cwd().name == "optimization" else None,
    ]

    for location in possible_locations:
        if location and location.exists() and (location / "src").exists():
            return location

    raise FileNotFoundError(
        "Could not find 'utag' or 'utag_paper' project root directory"
    )


try:
    main_dir = find_project_root()
    src_dir = main_dir / "src"

    if src_dir.exists():
        sys.path.insert(0, str(src_dir))
    else:
        raise FileNotFoundError(f"Source directory not found: {src_dir}")

    print(f"Project root: {main_dir}")
    print(f"Source directory: {src_dir}")
except Exception as e:
    print(f"Error setting up paths: {e}")
    sys.exit(1)

# try:
#     from imports import mi
# except ImportError as e:
#     print(f"Failed to import from imports.py: {e}")
#     sys.exit(1)
try:
    import microlive.microscopy as mi
except ImportError as e:
    print(f"Failed to import microlive.microscopy: {e}")
    sys.exit(1)

try:
    from tasep_models import (
        calculate_codon_elongation_rates,
        create_probe_vector,
        read_sequence,
        simulate_TASEP_SSA,
        tag_dict,
    )

    print("TASEP models loaded from tasep_models library")
except ImportError as e:
    print(f"Failed to import tasep_models: {e}")
    sys.exit(1)

# Input, output, and cache paths
ACF_DATA_FOLDER = main_dir / "optimization" / "source_files_ACF" / "data_ACF"
RESULTS_FOLDER = main_dir / "optimization" / "results_ACF"
CACHE_DIR = main_dir / "optimization" / "source_files_ACF" / "simulation_cache"

ACF_DATA_FOLDER.mkdir(parents=True, exist_ok=True)
RESULTS_FOLDER.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)


DATASET_CONFIGS = {
    "utag": {
        "tag_key": "U",
        "dna_filename": "pNZ208(pUB-24xUTagFullLength-KDM5B-MS2).dna",
    },
    "utag_c_free": {
        "tag_key": "U",
        "dna_filename": "pNZ208(pUB-24xUTagFullLength-KDM5B-MS2).dna",
    },
    "suntag": {
        "tag_key": "SUN",
        "dna_filename": "pNZ266(pUB-24xGCN4-KDM5B-MS2).dna",
    },
    "alfatag": {
        "tag_key": "ALFA",
        "dna_filename": "pNZ267 (pUB-24xALFAtag-KDM5B-MS2).dna",
    },
}


# Simulation cache
def get_simulation_hash(
    ki,
    global_elongation_rate,
    rna,
    first_probe_position_vector,
    burnin_time,
    gene_length,
    t_max,
    step_size_in_sec,
    number_repetitions,
    dataset_name,
):
    """Return a stable hash for one set of simulation parameters."""
    params = {
        "ki": float(ki),
        "global_elongation_rate": float(global_elongation_rate),
        "rna_sequence": str(rna),
        "first_probe_positions": list(first_probe_position_vector),
        "burnin_time": int(burnin_time),
        "gene_length": int(gene_length),
        "t_max": float(t_max),
        "step_size_in_sec": float(step_size_in_sec),
        "number_repetitions": int(number_repetitions),
        "dataset": str(dataset_name),
    }

    param_string = json.dumps(params, sort_keys=True)
    return hashlib.md5(param_string.encode()).hexdigest()


def save_simulation_cache(ssa_array, simulation_hash, cache_directory):
    """Save a simulated intensity array to the cache."""
    cache_file = cache_directory / f"simulation_{simulation_hash}.pkl"

    with open(cache_file, "wb") as file_handle:
        pickle.dump(ssa_array, file_handle)

    return cache_file


def load_simulation_cache(simulation_hash, cache_directory):
    """Load a simulated intensity array from the cache, if available."""
    cache_file = cache_directory / f"simulation_{simulation_hash}.pkl"

    if cache_file.exists():
        with open(cache_file, "rb") as file_handle:
            ssa_array = pickle.load(file_handle)
        return ssa_array

    return None


def print_cache_summary(cache_directory):
    """Print the number and total size of cached simulations."""
    cache_files = list(cache_directory.glob("simulation_*.pkl"))

    if cache_files:
        total_size_gb = sum(path.stat().st_size for path in cache_files) / 1024**3
        print(f"Cached simulations: {len(cache_files)}")
        print(f"Cache size: {total_size_gb:.2f} GB")
    else:
        print("No cached simulations found.")


def dataset_selection(dataset_name, downsample=False):
    """Load the experimental ACF and sequence configuration for a dataset."""
    if dataset_name not in DATASET_CONFIGS:
        raise ValueError(f"Unknown dataset: {dataset_name}")

    downsample_suffix = "_downsampled" if downsample else ""
    plot_name = f"{dataset_name}{downsample_suffix}"
    acf_data_file = ACF_DATA_FOLDER / f"df_ACF_{plot_name}.csv"

    config = DATASET_CONFIGS[dataset_name]
    tag_sequence = tag_dict[config["tag_key"]]
    dna_file_path = (
        main_dir / "gene_sequences" / "utag_project" / config["dna_filename"]
    )

    if not acf_data_file.exists():
        raise FileNotFoundError(f"ACF data file not found: {acf_data_file}")
    if not dna_file_path.exists():
        raise FileNotFoundError(f"DNA sequence file not found: {dna_file_path}")

    protein, rna, _, indexes_tags, _, _, _ = read_sequence(
        seq=dna_file_path, min_protein_length=50, TAG=[tag_sequence]
    )
    gene_length = len(protein) + 1  # Include the stop codon.
    tag_positions_first_probe_vector = indexes_tags[0]
    first_probe_position_vector = create_probe_vector(
        tag_positions_first_probe_vector, gene_length
    )

    acf_data = pd.read_csv(acf_data_file)
    data_mean_correlation = acf_data["mean_correlation"].values
    data_std_correlation = acf_data["std_correlation"].values
    data_lags = acf_data["lags"].values

    return (
        data_mean_correlation,
        data_std_correlation,
        data_lags,
        plot_name,
        rna,
        first_probe_position_vector,
        gene_length,
    )


def run_simulation(
    global_elongation_rate,
    ki,
    rna,
    first_probe_position_vector,
    burnin_time,
    gene_length,
    t_max,
    step_size_in_sec,
    number_repetitions,
    mad_threshold_factor,
    dataset_name,
    downsample=False,
    downsampling_factor=2,
    force_fresh=False,
):
    """Run a TASEP simulation and calculate its ACF."""
    sim_hash = get_simulation_hash(
        ki,
        global_elongation_rate,
        rna,
        first_probe_position_vector,
        burnin_time,
        gene_length,
        t_max,
        step_size_in_sec,
        number_repetitions,
        dataset_name,
    )

    ssa_array = None
    if not force_fresh:
        ssa_array = load_simulation_cache(sim_hash, CACHE_DIR)

    if ssa_array is None:
        codon_elongation_rates = calculate_codon_elongation_rates(
            rna, global_elongation_rate=global_elongation_rate
        )
        ssa_array = simulate_TASEP_SSA(
            ki,
            codon_elongation_rates,
            gene_length,
            t_max,
            time_interval_in_seconds=step_size_in_sec,
            number_repetitions=number_repetitions,
            first_probe_position_vector=first_probe_position_vector,
            burnin_time=burnin_time,
            constant_elongation_rate=global_elongation_rate,
            fast_output=True,
        )[2]

        save_simulation_cache(ssa_array, sim_hash, CACHE_DIR)

    (
        mean_correlation_ssa,
        std_correlation_ssa,
        lags_ssa,
        _,
        _,
    ) = mi.Correlation(
        primary_data=ssa_array,
        max_lag=None,
        nan_handling="forward_fill",
        shift_data=True,
        return_full=False,
        time_interval_between_frames_in_seconds=step_size_in_sec,
        use_bootstrap=False,
        show_plot=False,
        start_lag=0,
        fit_type="linear",
        de_correlation_threshold=0.01,
        correct_baseline=True,
        use_linear_projection_for_lag_0=True,
        save_plots=False,
        use_global_mean=False,
        remove_outliers=True,
        MAD_THRESHOLD_FACTOR=mad_threshold_factor,
        plot_individual_trajectories=False,
        y_axes_min_max_list_values=None,
        x_axes_min_max_list_values=None,
        plot_title=None,
    ).run()

    if downsample:
        mean_correlation_ssa = mean_correlation_ssa[::downsampling_factor]
        std_correlation_ssa = std_correlation_ssa[::downsampling_factor]
        lags_ssa = lags_ssa[::downsampling_factor]

    return mean_correlation_ssa, std_correlation_ssa, lags_ssa


def comparing_acfs(
    mean_correlation,
    mean_correlation_ssa,
    err_correlation,
    start_lag=1,
    end_lag=150,
):
    """Return the chi-squared-like difference between two ACFs."""
    err_correlation = err_correlation + 1e-6
    end_lag = min(end_lag, len(mean_correlation), len(mean_correlation_ssa))

    residuals = (
        mean_correlation[start_lag:end_lag] - mean_correlation_ssa[start_lag:end_lag]
    )
    objective_value = np.sum(residuals**2 / err_correlation[start_lag:end_lag] ** 2)
    return objective_value


def optimization_function(parameters):
    """Evaluate one pair of initiation and elongation rates."""
    ki, global_elongation_rate = parameters

    mean_correlation_ssa, _, _ = run_simulation(
        global_elongation_rate=global_elongation_rate,
        ki=ki,
        rna=rna,
        first_probe_position_vector=first_probe_position_vector,
        burnin_time=BURN_IN_TIME_SEC,
        gene_length=gene_length,
        t_max=SIMULATION_DURATION_SEC,
        step_size_in_sec=STEP_SIZE_IN_SEC,
        number_repetitions=NUMBER_REPETITIONS,
        mad_threshold_factor=MAD_THRESHOLD_FACTOR,
        dataset_name=dataset,
        downsample=False,
        downsampling_factor=DOWNSAMPLING_FACTOR,
        force_fresh=False,
    )

    return comparing_acfs(
        mean_correlation,
        mean_correlation_ssa,
        std_correlation,
        start_lag=1,
        end_lag=130,
    )


def plot_correlation_comparison(
    mean_correlation,
    std_correlation,
    lags,
    mean_correlation_ssa=None,
    std_correlation_ssa=None,
    lags_ssa=None,
    plot_name="temp",
    start_lag=0,
    y_axes_min_max_list_values=None,
    x_axes_min_max_list_values=None,
    show_legend=True,
):
    """Plot the experimental and simulated ACFs."""
    fig, ax = plt.subplots(figsize=(5, 3), facecolor="white")

    lags_minutes = lags / 60
    max_lag = len(lags)

    ax.plot(
        lags_minutes[start_lag:max_lag],
        mean_correlation[start_lag:max_lag],
        "o-",
        color="dimgray",
        linewidth=1,
        label="Experimental",
        alpha=0.5,
    )
    ax.fill_between(
        lags_minutes[start_lag:max_lag],
        mean_correlation[start_lag:max_lag] - std_correlation[start_lag:max_lag],
        mean_correlation[start_lag:max_lag] + std_correlation[start_lag:max_lag],
        color="dimgray",
        alpha=0.5,
    )

    if mean_correlation_ssa is not None:
        lags_ssa_minutes = lags_ssa / 60
        ax.plot(
            lags_ssa_minutes[start_lag:max_lag],
            mean_correlation_ssa[start_lag:max_lag],
            "o-",
            color="goldenrod",
            linewidth=1,
            label="Simulation",
            alpha=0.5,
        )
        ax.fill_between(
            lags_ssa_minutes[start_lag:max_lag],
            mean_correlation_ssa[start_lag:max_lag]
            - std_correlation_ssa[start_lag:max_lag],
            mean_correlation_ssa[start_lag:max_lag]
            + std_correlation_ssa[start_lag:max_lag],
            color="goldenrod",
            alpha=0.5,
        )

    if show_legend:
        ax.legend(loc="upper right")

    if x_axes_min_max_list_values is not None:
        ax.set_xlim(x_axes_min_max_list_values)
    if y_axes_min_max_list_values is not None:
        ax.set_ylim(y_axes_min_max_list_values)

    ax.set_xlabel(r"$\tau$ (min)", fontsize=14)
    ax.set_ylabel(r"$G(\tau)$", fontsize=14)

    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1.5)

    output_file = RESULTS_FOLDER / f"acf_comp_{plot_name}.png"
    fig.tight_layout()
    fig.savefig(output_file, dpi=600)
    plt.close(fig)

    print(f"Correlation comparison plot saved: {output_file}")


def grid_search_ki_ke(
    optimization_function,
    dim_size=10,
    ki_min=0.01,
    ki_max=1.0,
    ke_min=1.0,
    ke_max=10.0,
    n_jobs=-1,
    top_percentile=5.0,
    confidence_level=0.95,
):
    """Evaluate a regular ki-ke grid and summarize its best solutions."""
    ki_values = np.linspace(ki_min, ki_max, dim_size)
    ke_values = np.linspace(ke_min, ke_max, dim_size)
    param_grid = [(float(ki), float(ke)) for ki in ki_values for ke in ke_values]

    def compute_cost(ki_ke):
        ki_value, ke_value = ki_ke
        try:
            cost = optimization_function([ki_value, ke_value])
            if not np.isfinite(cost):
                return ki_value, ke_value, 1e10
            return ki_value, ke_value, cost
        except Exception:
            return ki_value, ke_value, 1e10

    print(f"Running {len(param_grid)} simulations with {n_jobs} workers.")
    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_cost)(pair) for pair in param_grid
    )

    results_df = pd.DataFrame(results, columns=["ki", "ke", "cost"])
    valid_mask = np.isfinite(results_df["cost"]) & (results_df["cost"] < 1e9)
    results_df_clean = results_df[valid_mask].copy()

    valid_percentage = 100 * len(results_df_clean) / len(results_df)
    failed_count = len(results_df) - len(results_df_clean)
    print(
        f"Valid results: {len(results_df_clean)}/{len(results_df)} "
        f"({valid_percentage:.1f}%)"
    )
    if failed_count:
        print(f"Failed evaluations: {failed_count}")

    if len(results_df_clean) == 0:
        raise ValueError("No valid optimization results found!")

    best_idx = results_df_clean["cost"].idxmin()
    best_cost = results_df_clean.loc[best_idx, "cost"]

    n_top = max(1, int(len(results_df_clean) * top_percentile / 100))
    top_solutions = results_df_clean.nsmallest(n_top, "cost")

    robust_ki = np.median(top_solutions["ki"])
    robust_ke = np.median(top_solutions["ke"])

    alpha = 1 - confidence_level
    low_percentile = 100 * alpha / 2
    high_percentile = 100 * (1 - alpha / 2)

    ki_ci = np.percentile(top_solutions["ki"], [low_percentile, high_percentile])
    ke_ci = np.percentile(top_solutions["ke"], [low_percentile, high_percentile])

    return results_df_clean, best_cost, robust_ki, robust_ke, ki_ci, ke_ci


def print_section(title):
    """Print a consistent heading in the terminal log."""
    separator = "=" * 70
    print(f"\n{separator}\n{title}\n{separator}")


def save_optimization_results(
    results_df,
    best_cost,
    robust_ki,
    robust_ke,
    ki_ci,
    ke_ci,
    plot_name,
    dataset_name,
):
    """Save the grid-search results and parameter summary."""
    parameters = {
        "ki": robust_ki,
        "ke": robust_ke,
        "ki_90%_CI_min": ki_ci[0],
        "ki_90%_CI_max": ki_ci[1],
        "ke_90%_CI_min": ke_ci[0],
        "ke_90%_CI_max": ke_ci[1],
        "best_cost": best_cost,
    }
    df_parameters = pd.DataFrame(parameters.items(), columns=["parameter", "value"])
    param_file = RESULTS_FOLDER / f"optimized_params_{plot_name}.csv"
    df_parameters.to_csv(param_file, index=False)

    results_file = RESULTS_FOLDER / f"results_df_{dataset_name}.csv"
    results_df.to_csv(results_file, index=False)

    print(f"\nOptimized parameters saved to: {param_file}")
    print(f"Full grid-search results saved to: {results_file}")


def main():
    """Run the complete optimization workflow for the selected dataset."""
    global mean_correlation
    global std_correlation
    global rna
    global first_probe_position_vector
    global gene_length

    print_section("TASEP OPTIMIZATION")
    print(f"Selected dataset: {dataset} (index {dataset_index})")
    print(f"Parallel workers: {PARALLEL_JOBS}")
    print(f"Cache directory: {CACHE_DIR}")
    print_cache_summary(CACHE_DIR)

    print(f"\nLoading dataset: {dataset}")
    (
        mean_correlation,
        std_correlation,
        lags,
        plot_name,
        rna,
        first_probe_position_vector,
        gene_length,
    ) = dataset_selection(dataset_name=dataset, downsample=DOWNSAMPLE)

    print_section("GRID SEARCH")
    print(f"Grid: {GRID_SIZE} x {GRID_SIZE} = {GRID_SIZE**2} evaluations")
    print(f"ki range: [{KI_MIN}, {KI_MAX}] s^-1")
    print(f"ke range: [{KE_MIN}, {KE_MAX}] aa/s")

    results_df, best_cost, robust_ki, robust_ke, ki_ci, ke_ci = grid_search_ki_ke(
        optimization_function=optimization_function,
        dim_size=GRID_SIZE,
        ki_min=KI_MIN,
        ki_max=KI_MAX,
        ke_min=KE_MIN,
        ke_max=KE_MAX,
        n_jobs=PARALLEL_JOBS,
        top_percentile=TOP_PERCENTILE,
        confidence_level=CONFIDENCE_LEVEL,
    )

    print_section("OPTIMIZATION RESULTS")
    print(f"Robust estimates (median of top {TOP_PERCENTILE}%):")
    print(f"  ki = {robust_ki:.6f} s^-1; " f"90% CI: [{ki_ci[0]:.6f}, {ki_ci[1]:.6f}]")
    print(f"  ke = {robust_ke:.6f} aa/s; " f"90% CI: [{ke_ci[0]:.6f}, {ke_ci[1]:.6f}]")
    print(f"  Minimum grid cost: {best_cost:.8f}")

    save_optimization_results(
        results_df=results_df,
        best_cost=best_cost,
        robust_ki=robust_ki,
        robust_ke=robust_ke,
        ki_ci=ki_ci,
        ke_ci=ke_ci,
        plot_name=plot_name,
        dataset_name=dataset,
    )

    print_section("FINAL SIMULATION WITH OPTIMIZED PARAMETERS")
    (
        mean_correlation_ssa_final,
        std_correlation_ssa_final,
        lags_ssa_final,
    ) = run_simulation(
        global_elongation_rate=robust_ke,
        ki=robust_ki,
        rna=rna,
        first_probe_position_vector=first_probe_position_vector,
        burnin_time=BURN_IN_TIME_SEC,
        gene_length=gene_length,
        t_max=SIMULATION_DURATION_SEC,
        step_size_in_sec=STEP_SIZE_IN_SEC,
        number_repetitions=NUMBER_REPETITIONS,
        mad_threshold_factor=MAD_THRESHOLD_FACTOR,
        dataset_name=dataset,
        downsample=DOWNSAMPLE,
        downsampling_factor=DOWNSAMPLING_FACTOR,
        force_fresh=True,
    )

    if mean_correlation_ssa_final is None or len(mean_correlation_ssa_final) == 0:
        print("ERROR: The final simulation did not produce results.")
    else:
        print("Final simulation completed successfully.")
        print(f"  Result shape: {mean_correlation_ssa_final.shape}")
        print(
            "  Mean correlation range: "
            f"[{np.min(mean_correlation_ssa_final):.4f}, "
            f"{np.max(mean_correlation_ssa_final):.4f}]"
        )

    print("\nGenerating the final comparison plot.")
    plot_correlation_comparison(
        mean_correlation=mean_correlation,
        std_correlation=std_correlation,
        lags=lags,
        mean_correlation_ssa=mean_correlation_ssa_final,
        std_correlation_ssa=std_correlation_ssa_final,
        lags_ssa=lags_ssa_final,
        plot_name=f"{plot_name}_final",
        start_lag=0,
        y_axes_min_max_list_values=Y_AXIS_LIMITS,
        x_axes_min_max_list_values=X_AXIS_LIMITS,
        show_legend=True,
    )

    print_section("EXECUTION COMPLETE")
    print_cache_summary(CACHE_DIR)
    print(f"All results saved to: {RESULTS_FOLDER}\n")


if __name__ == "__main__":
    main()
