#!/usr/bin/env python3
import re
from collections.abc import Mapping
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

current_dir = Path().resolve()
from microlive.imports import *
from microlive import microscopy as mi
import tasep_models as tm
from tasep_models import *


plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'axes.edgecolor': 'black',
        'axes.linewidth': 1.5,
        'font.family': 'Arial',
        'font.sans-serif': ['Arial'],
        'axes.labelsize': 16,
        'axes.titlesize': 16,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'axes.labelcolor': 'black',
        'text.color': 'black',
        'xtick.color': 'black',
        'ytick.color': 'black',
    })

def calculate_number_of_particles_per_frame(particle_counts_per_frame, inhibitor_frame_index):
    """Normalize the number of particles per frame to the average before treatment.

    All frames are divided by the average particle count before treatment.
    If that average is zero, returns an array of zeros.

    Args:
        particle_counts_per_frame: 1D array of particle counts per frame.
        inhibitor_frame_index: The frame index at which treatment starts.

    Returns:
        Tuple of (normalized_particles, average_particles_before_treatment).
    """
    # Compute the average particle count before treatment
    pre_counts = particle_counts_per_frame[:inhibitor_frame_index]
    average_particles_before_treatment = pre_counts.mean()

    # Normalize all frames by the pre-treatment average
    if average_particles_before_treatment == 0:
        normalized_particles = np.zeros_like(particle_counts_per_frame, dtype=float)
    else:
        normalized_particles = particle_counts_per_frame / average_particles_before_treatment

    return normalized_particles, average_particles_before_treatment


def _resolve_experiment_lif_files(experiment_lif_files, experiment_labels=None):
    """Resolve labeled LIF inputs from a mapping, file list, or directory."""
    if isinstance(experiment_lif_files, Mapping):
        labeled_paths = [
            (str(label), Path(file_path))
            for label, file_path in experiment_lif_files.items()
        ]
    else:
        if isinstance(experiment_lif_files, (str, Path)):
            input_path = Path(experiment_lif_files)
            lif_paths = sorted(input_path.glob('*.lif')) if input_path.is_dir() else [input_path]
        else:
            lif_paths = [Path(file_path) for file_path in experiment_lif_files]

        if experiment_labels is None:
            labeled_paths = [(file_path.stem, file_path) for file_path in lif_paths]
        elif isinstance(experiment_labels, Mapping):
            labeled_paths = []
            for file_path in lif_paths:
                label = experiment_labels.get(
                    file_path.stem,
                    experiment_labels.get(file_path.name, str(file_path.stem)),
                )
                labeled_paths.append((str(label), file_path))
        else:
            if len(experiment_labels) != len(lif_paths):
                raise ValueError('experiment_labels must match the number of LIF files.')
            labeled_paths = list(zip(map(str, experiment_labels), lif_paths))

    if not labeled_paths:
        raise ValueError('No LIF files were provided or discovered.')

    resolved_paths = {}
    stems = set()
    for label, file_path in labeled_paths:
        file_path = file_path.expanduser().resolve()
        if not file_path.is_file() or file_path.suffix.lower() != '.lif':
            raise FileNotFoundError(f'Invalid LIF file for {label!r}: {file_path}')
        if label in resolved_paths:
            raise ValueError(f'Duplicate experiment label: {label!r}')
        if file_path.stem in stems:
            raise ValueError(f'Duplicate LIF stem: {file_path.stem!r}')
        resolved_paths[label] = file_path
        stems.add(file_path.stem)
    return resolved_paths


def _resolve_results_dirs(data_dir):
    """Resolve one or more directories containing result subfolders."""
    if isinstance(data_dir, Mapping):
        directory_inputs = list(data_dir.values())
    elif isinstance(data_dir, (str, Path)):
        directory_inputs = [data_dir]
    else:
        directory_inputs = list(data_dir)
    if not directory_inputs:
        raise ValueError('At least one results directory is required.')

    results_dirs = []
    for directory_input in directory_inputs:
        results_dir = Path(directory_input).expanduser().resolve()
        if not results_dir.is_dir():
            raise NotADirectoryError(f'Results directory not found: {results_dir}')
        if results_dir not in results_dirs:
            results_dirs.append(results_dir)
    return results_dirs


def _read_microlive_metadata(metadata_path):
    """Read image-specific timing and source information from MicroLive metadata."""
    metadata_text = Path(metadata_path).read_text(encoding='utf-8', errors='replace')

    def get_value(field_name):
        pattern = rf'^\s*{re.escape(field_name)}\.+\s*(.*?)\s*$'
        match = re.search(pattern, metadata_text, flags=re.MULTILINE)
        if match is None:
            raise ValueError(f'{field_name!r} was not found in {metadata_path}.')
        return match.group(1)

    frame_interval_sec = float(get_value('Time Interval (s)'))
    active_frame_count = int(get_value('Active Frame Count'))
    if frame_interval_sec <= 0 or active_frame_count <= 0:
        raise ValueError(f'Invalid timing metadata in {metadata_path}.')
    return {
        'source_lif_path': Path(get_value('Data Folder Path')),
        'selected_image': get_value('Image Name'),
        'frame_interval_sec': frame_interval_sec,
        'active_frame_count': active_frame_count,
    }


def load_particle_counts_by_experiment(
        data_dir, experiment_lif_files, experiment_labels=None,
        condition_labels=None, baseline_frames=1, particle_column='particle',
        after_only=True, verbose=False):
    """Load and normalize particle counts for labeled microscopy experiments.

    Args:
        data_dir: Directory containing ``results_*`` subfolders, or a
            sequence/mapping of multiple results directories.
        experiment_lif_files: Mapping of experiment label to LIF path, a list
            of LIF paths, a single LIF path, or a directory containing LIFs.
        experiment_labels: Optional labels for non-mapping LIF inputs. May be
            a stem-to-label mapping or a sequence matching the LIF file order.
        condition_labels: Optional experiment-to-display-label mapping.
        baseline_frames: Number of initial frames used for normalization.
            The default of 1 normalizes every cell to its first frame.
        particle_column: Particle identifier used for per-frame unique counts.
        after_only: If True, exclude images explicitly named ``Before ActD``.
        verbose: If True, print an experiment summary.

    Returns:
        Tidy DataFrame with raw and normalized counts for every cell and frame.
    """
    data_dirs = _resolve_results_dirs(data_dir)
    if not isinstance(baseline_frames, int) or baseline_frames < 1:
        raise ValueError('baseline_frames must be a positive integer.')

    lif_files = _resolve_experiment_lif_files(
        experiment_lif_files, experiment_labels=experiment_labels,
    )
    experiment_by_lif_stem = {
        file_path.stem: experiment for experiment, file_path in lif_files.items()
    }
    condition_labels = condition_labels or {}
    result_folders = sorted(
        (results_dir, folder)
        for results_dir in data_dirs
        for folder in results_dir.iterdir()
        if folder.is_dir() and folder.name.startswith('results_')
    )
    if not result_folders:
        joined_dirs = ', '.join(map(str, data_dirs))
        raise ValueError(f'No results_* folders found in: {joined_dirs}')

    dataframe_parts = []
    skipped_before_images = []
    unmatched_folders = []
    for results_dir, result_folder in result_folders:
        tracking_files = sorted(result_folder.glob('tracking_*.csv'))
        metadata_files = sorted(result_folder.glob('Metadata_*.txt'))
        if len(tracking_files) != 1 or len(metadata_files) != 1:
            raise ValueError(
                f'{result_folder.name} must contain exactly one tracking_*.csv '
                f'and one Metadata_*.txt file; found {len(tracking_files)} and '
                f'{len(metadata_files)}, respectively.'
            )

        metadata = _read_microlive_metadata(metadata_files[0])
        source_lif_stem = metadata['source_lif_path'].stem
        experiment = experiment_by_lif_stem.get(source_lif_stem)
        if experiment is None:
            unmatched_folders.append(str(result_folder))
            continue
        selected_image_lower = metadata['selected_image'].lower().replace('_', ' ')
        if after_only and 'before actd' in selected_image_lower:
            skipped_before_images.append(result_folder.name)
            continue
        if baseline_frames > metadata['active_frame_count']:
            raise ValueError(
                f'baseline_frames={baseline_frames} exceeds the active frame '
                f'count in {metadata_files[0].name}.'
            )

        tracking_df = pd.read_csv(tracking_files[0])
        required_columns = {'cell_id', 'frame'}
        missing_columns = required_columns.difference(tracking_df.columns)
        if missing_columns:
            raise ValueError(
                f'{tracking_files[0].name} is missing columns: '
                f'{sorted(missing_columns)}'
            )
        count_column = particle_column
        if count_column not in tracking_df.columns:
            count_column = 'unique_particle' if 'unique_particle' in tracking_df.columns else None
        if count_column is None:
            raise ValueError(
                f'{tracking_files[0].name} has neither {particle_column!r} nor '
                f'"unique_particle".'
            )
        if tracking_df.empty or tracking_df['cell_id'].isna().any():
            raise ValueError(f'{tracking_files[0].name} contains no valid cells.')

        tracking_df = tracking_df.copy()
        numeric_frames = pd.to_numeric(tracking_df['frame'], errors='raise')
        if not np.allclose(numeric_frames, np.round(numeric_frames)):
            raise ValueError(f'{tracking_files[0].name} contains non-integer frames.')
        tracking_df['frame'] = numeric_frames.astype(int)
        active_frame_count = metadata['active_frame_count']
        if ((tracking_df['frame'] < 0) | (tracking_df['frame'] >= active_frame_count)).any():
            raise ValueError(
                f'{tracking_files[0].name} contains frames outside its metadata range.'
            )

        frame_interval_sec = metadata['frame_interval_sec']
        if 'time' in tracking_df.columns:
            observed_times = pd.to_numeric(tracking_df['time'], errors='coerce')
            expected_times = tracking_df['frame'].to_numpy() * frame_interval_sec
            valid_times = observed_times.notna().to_numpy()
            if valid_times.any() and not np.allclose(
                    observed_times.to_numpy()[valid_times], expected_times[valid_times],
                    rtol=0, atol=max(1e-6, frame_interval_sec * 1e-6)):
                raise ValueError(
                    f'Tracking times disagree with selected-image metadata in '
                    f'{tracking_files[0].name}.'
                )

        frame_index = pd.Index(range(active_frame_count), name='frame')
        for cell_id, cell_df in tracking_df.groupby('cell_id', sort=True):
            raw_counts = (
                cell_df.groupby('frame')[count_column]
                .nunique(dropna=True)
                .reindex(frame_index, fill_value=0)
                .to_numpy(dtype=float)
            )
            normalized_counts, baseline_count = calculate_number_of_particles_per_frame(
                raw_counts, baseline_frames,
            )
            cell_key = f'{experiment}::{result_folder.name}::cell_{cell_id}'
            dataframe_parts.append(pd.DataFrame({
                'experiment': experiment,
                'condition': condition_labels.get(experiment, experiment),
                'source_lif': str(lif_files[experiment]),
                'results_dir': str(results_dir),
                'result_folder': result_folder.name,
                'selected_image': metadata['selected_image'],
                'cell_id': cell_id,
                'cell_key': cell_key,
                'frame': frame_index.to_numpy(),
                'time_sec': frame_index.to_numpy(dtype=float) * frame_interval_sec,
                'time_min': frame_index.to_numpy(dtype=float) * frame_interval_sec / 60.0,
                'frame_interval_sec': frame_interval_sec,
                'active_frame_count': active_frame_count,
                'particle_count': raw_counts,
                'normalized_particle_count': normalized_counts,
                'baseline_particle_count': baseline_count,
            }))

    if unmatched_folders:
        names = '\n  '.join(unmatched_folders)
        raise ValueError(f'Result folders did not match a supplied LIF:\n  {names}')
    if not dataframe_parts:
        raise ValueError('No after-treatment/control tracking data were loaded.')

    particle_counts_df = pd.concat(dataframe_parts, ignore_index=True)
    experiment_order = {name: index for index, name in enumerate(lif_files)}
    particle_counts_df['_experiment_order'] = particle_counts_df['experiment'].map(experiment_order)
    particle_counts_df = (
        particle_counts_df
        .sort_values(['_experiment_order', 'cell_key', 'frame'])
        .drop(columns='_experiment_order')
        .reset_index(drop=True)
    )
    particle_counts_df.attrs['experiment_order'] = list(lif_files)
    particle_counts_df.attrs['results_dirs'] = list(map(str, data_dirs))
    particle_counts_df.attrs['skipped_before_images'] = skipped_before_images

    if verbose:
        summary = particle_counts_df.groupby('experiment', sort=False).agg(
            cells=('cell_key', 'nunique'),
            frames=('active_frame_count', 'max'),
            interval_sec=('frame_interval_sec', 'first'),
            max_time_min=('time_min', 'max'),
        )
        print(summary.to_string())
        if skipped_before_images:
            print(f'Skipped {len(skipped_before_images)} Before ActD result folders.')
    return particle_counts_df


def calculate_intensity(particle_counts_per_frame, sum_intensities_per_frame, inhibitor_frame_index, normalization_method='mean', percentile_range=(5, 95)):
    """Normalize the intensity per frame.

    For frames before the treatment, each frame's intensity is given by
    sum_intensities / particle_counts. If the particle count is zero in a frame,
    the normalized intensity is set to zero.

    For frames after the treatment, the sum intensities are divided by the average
    particle count before treatment. If that average is zero, zeros are returned for
    all frames after treatment.

    Args:
        particle_counts_per_frame: 1D array of particle counts per frame.
        sum_intensities_per_frame: 1D array of sum intensities per frame.
        inhibitor_frame_index: The frame index at which treatment starts.
        normalization_method: 'mean' (default) divides by pre-treatment mean,
            'minmax' scales the trajectory to [0, 1] range,
            'percentile' scales using percentile bounds (robust to outliers).
            None: no normalization.
        percentile_range: Percentile bounds for 'percentile' method.
            Default (5, 95). Use (1, 99) for wider range.

    Returns:
        Tuple of (normalized_intensities, raw_avg_intensities, avg_particles_before_treatment).
    """
    # Compute the average particle count before treatment.
    pre_counts = particle_counts_per_frame[:inhibitor_frame_index]
    average_particles_before_treatment = pre_counts.mean()

    # For frames before treatment, avoid division by zero:
    pre_intensities = sum_intensities_per_frame[:inhibitor_frame_index]
    intensity_before_treatment = np.divide(
        pre_intensities, pre_counts,
        out=np.zeros(pre_intensities.shape, dtype=float),
        where=pre_counts != 0,
    )
    # For frames after treatment, if the average is zero then return zeros.
    post_intensities = sum_intensities_per_frame[inhibitor_frame_index:]
    if average_particles_before_treatment == 0:
        intensity_after_treatment = np.zeros_like(post_intensities)
    else:
        intensity_after_treatment = post_intensities / average_particles_before_treatment
    # Combine the two segments and return the result.
    average_intensity_with_respect_number_particles = np.concatenate([intensity_before_treatment, intensity_after_treatment])

    # Apply normalization
    # 'minmax' and 'percentile' are handled globally in process_inhibitor_data
    # (across all cells), so here they pass through raw per-cell intensities.
    # 'mean' is inherently per-cell (divides by pre-treatment mean).
    if normalization_method == 'mean':
        mean_before_treatment = average_intensity_with_respect_number_particles[:inhibitor_frame_index].mean()
        if mean_before_treatment == 0:
            intensities_normalized_before_treatment_intensity = np.zeros_like(average_intensity_with_respect_number_particles)
        else:
            intensities_normalized_before_treatment_intensity = average_intensity_with_respect_number_particles / mean_before_treatment
    else:
        intensities_normalized_before_treatment_intensity = average_intensity_with_respect_number_particles

    return intensities_normalized_before_treatment_intensity, average_intensity_with_respect_number_particles, average_particles_before_treatment


# ── Model definitions ────────────────────────────────────────────────────────

def _linear_model(x, a, b):
    """Linear decay: y = a*x + b"""
    return a * x + b


def _exponential_model(x, A, tau, C):
    """Exponential decay: y = A * exp(-x/τ) + C"""
    return A * np.exp(-x / tau) + C


def _heaviside_model(x, A, T, C):
    """Larson 2011 Heaviside-ramp: y = A*(1 - x/T)*H(T - x) + C

    Linear decay from A+C to C, then flat at C for x > T.
    """
    x = np.asarray(x, dtype=float)
    return np.where(x <= T, A * (1.0 - x / T) + C, C)


# ── Fitting function ─────────────────────────────────────────────────────────

def fit_inhibitor_model(x_data, y_data, err_data=None, model='exponential',
                        fit_start_idx=None, fit_end_idx=None,
                        runoff_fraction=0.95, basal_value=None):
    """Fit inhibitor run-off data to a decay model.

    Args:
        x_data: Time array (e.g., time in minutes, recentered so 0 = inhibitor).
        y_data: Mean intensity trajectory (1D).
        err_data: Per-point measurement uncertainty (e.g., SEM from individual
            cells). Same length as y_data. When provided, curve_fit performs
            weighted least-squares and χ² is computed as Σ[(y-f)²/σ²].
            When None, unweighted fitting is used.
        model: One of 'linear', 'exponential', 'heaviside',
            'linear_extrapolated'.
        basal_value: Baseline intensity for 'linear_extrapolated' model.
            When provided, t_runoff is the time at which the fitted line
            crosses this value. Typically computed from the last
            background_frames of the normalized mean trajectory. When None
            and model is 'linear_extrapolated', falls back to estimating
            baseline from the last 20% of the full y_data array.
        fit_start_idx: Index into x_data/y_data for the start of the fitting
            range. Defaults to 0 (start of the array).
        fit_end_idx: Index into x_data/y_data for the end of the fitting range
            (inclusive). Defaults to len(x_data) - 1 (end of the array).
        runoff_fraction: Fraction of total decay used to define run-off time
            (default 0.95).

    Returns:
        On success, a dictionary with keys: 'model', 'params',
        'fitted_curve', 't_half', 't_runoff', 'R2', 'chi2',
        'chi2_reduced', 'dof'. Returns None if fitting fails.
    """
    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)
    if err_data is not None:
        err_data = np.asarray(err_data, dtype=float)

    # Default range: full array
    i0 = fit_start_idx if fit_start_idx is not None else 0
    i1 = (fit_end_idx + 1) if fit_end_idx is not None else len(x_data)

    x_fit = x_data[i0:i1]
    y_fit = y_data[i0:i1]
    err_fit = err_data[i0:i1] if err_data is not None else None

    # Remove NaN values (e.g., from artifact removal at inhibitor frame)
    valid = np.isfinite(x_fit) & np.isfinite(y_fit)
    if err_fit is not None:
        valid = valid & np.isfinite(err_fit)
    x_fit = x_fit[valid]
    y_fit = y_fit[valid]
    if err_fit is not None:
        err_fit = err_fit[valid]
        # Replace zero uncertainties with a small fraction of the minimum
        # non-zero value. When ALL uncertainties are zero (identical replicates
        # at every time point — cannot happen with real multi-cell data), this
        # falls back to a fixed 1e-10 so curve_fit still converges.
        err_fit = np.where(err_fit == 0, np.min(err_fit[err_fit > 0]) * 0.1 if np.any(err_fit > 0) else 1e-10, err_fit)
        sigma_kwarg = {'sigma': err_fit, 'absolute_sigma': True}
    else:
        sigma_kwarg = {}

    if len(x_fit) < 3:
        print('fit_inhibitor_model: not enough data points to fit.')
        return None

    model = model.lower().strip()

    try:
        if model == 'linear':
            # y = a*x + b
            popt, pcov = curve_fit(_linear_model, x_fit, y_fit, **sigma_kwarg)
            a, b = popt
            perr = np.sqrt(np.diag(pcov))
            fitted_full = _linear_model(x_data, *popt)

            # Derived quantities
            # Estimate actual baseline from last 20% of data
            tail = max(1, len(y_fit) // 5)
            Iss = float(np.mean(y_fit[-tail:]))
            I0 = b  # intensity at x = 0 (fitted intercept)
            if a != 0 and I0 != Iss:
                t_half = (I0 - (I0 + Iss) / 2.0) / (-a)   # when y = midpoint
                t_runoff = 2.0 * t_half
            else:
                t_half = np.inf
                t_runoff = np.inf

            params = {'a (slope)': a, 'b (intercept)': b,
                      'Iss (baseline)': Iss,
                      'a_err': perr[0], 'b_err': perr[1]}

        elif model == 'exponential':
            # y = A * exp(-x/τ) + C
            # Initial guesses (robust to negative/positive baselines)
            tail = max(1, len(y_fit) // 5)
            C0 = max(float(np.mean(y_fit[-tail:])), 0.0)  # nonnegative baseline
            A0 = max(float(y_fit[0]) - C0, 1e-8)
            tau0 = (x_fit[-1] - x_fit[0]) / 3.0
            popt, pcov = curve_fit(
                _exponential_model, x_fit, y_fit,
                p0=[A0, tau0, C0],
                # Normalized spot counts cannot have a negative plateau.
                # Constraining C prevents a poorly identified slow decay
                # from producing negative half-levels and enormous half-lives.
                bounds=([0, 1e-6, 0], [np.inf, np.inf, np.inf]),
                maxfev=50000,
                **sigma_kwarg,
            )
            A, tau, C = popt
            perr = np.sqrt(np.diag(pcov))
            fitted_full = _exponential_model(x_data, *popt)

            # Derived quantities
            t_half = tau * np.log(2)
            t_runoff = 2.0 * t_half

            params = {'A (amplitude)': A, 'tau (time constant)': tau,
                      'C (baseline)': C,
                      'A_err': perr[0], 'tau_err': perr[1], 'C_err': perr[2]}

        elif model == 'heaviside':
            # y = A * (1 - x/T) * H(T - x) + C   (Larson 2011)
            # Initial guesses (robust to negative/positive baselines)
            tail = max(1, len(y_fit) // 5)
            C0 = float(np.mean(y_fit[-tail:]))           # baseline from last 20%
            A0 = max(float(y_fit[0]) - C0, 1e-8)         # amplitude above baseline
            # Smart T guess: find where data first drops to baseline level
            crossings = np.where(y_fit <= C0)[0]
            if len(crossings) > 0:
                T0 = float(x_fit[crossings[0]] - x_fit[0])
            else:
                T0 = float((x_fit[-1] - x_fit[0]) / 2.0)  # fallback: half the range
            T0 = max(T0, 1.0)  # at least 1 minute
            popt, pcov = curve_fit(
                _heaviside_model, x_fit, y_fit,
                p0=[A0, T0, C0],
                bounds=([0, 1e-6, -np.inf], [np.inf, np.inf, np.inf]),
                maxfev=50000,
                **sigma_kwarg,
            )
            A, T, C = popt
            perr = np.sqrt(np.diag(pcov))
            fitted_full = _heaviside_model(x_data, *popt)

            # Derived quantities
            t_half = T / 2.0
            t_runoff = 2.0 * t_half  # = T (the full dwell time)

            params = {'A (amplitude)': A, 'T (dwell/run-off time)': T,
                      'C (baseline)': C,
                      'A_err': perr[0], 'T_err': perr[1], 'C_err': perr[2]}

        elif model == 'linear_extrapolated':
            # y = a*x + b  (same linear fit, but t_runoff = basal crossing)
            popt, pcov = curve_fit(_linear_model, x_fit, y_fit, **sigma_kwarg)
            a, b = popt
            perr = np.sqrt(np.diag(pcov))
            fitted_full = _linear_model(x_data, *popt)

            # Use basal value from caller; fallback to tail of full data
            if basal_value is not None:
                Iss = float(basal_value)
            else:
                tail = max(1, len(y_data) // 5)
                Iss = float(np.mean(y_data[-tail:]))

            I0 = b  # intensity at x = 0 (fitted intercept)

            # t_runoff = time where fit line crosses basal
            # Solve: a * t + b = Iss  →  t = (Iss - b) / a
            if a != 0 and I0 != Iss:
                t_runoff = (Iss - b) / a
                t_half = t_runoff / 2.0
            else:
                t_half = np.inf
                t_runoff = np.inf

            params = {'a (slope)': a, 'b (intercept)': b,
                      'Iss (baseline)': Iss,
                      'a_err': perr[0], 'b_err': perr[1]}

        else:
            print(f'fit_inhibitor_model: unknown model "{model}". '
                  f'Choose from: linear, exponential, heaviside, '
                  f'linear_extrapolated.')
            return None

        # Goodness-of-fit metrics (computed on the fitting window only)
        n_params = len(popt)
        n_data = len(y_fit)
        dof = n_data - n_params
        y_pred_fit = fitted_full[i0:i1][valid]
        residuals = y_fit - y_pred_fit

        # Unweighted sums of squares (always computed)
        ss_res = float(np.sum(residuals ** 2))
        ss_tot = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot != 0 else np.nan

        # Chi-squared: weighted if err_data provided, unweighted otherwise
        if err_fit is not None:
            chi2 = float(np.sum((residuals / err_fit) ** 2))
        else:
            chi2 = ss_res  # equivalent to unweighted χ²
        chi2_red = chi2 / dof if dof > 0 else np.nan

        result = {
            'model': model,
            'params': params,
            'fitted_curve': fitted_full,
            't_half': t_half,
            't_runoff': t_runoff,
            'runoff_fraction': runoff_fraction,
            'chi2': chi2,
            'chi2_reduced': chi2_red,
            'dof': dof,
            'R2': r_squared,
            'n_data': n_data,
        }
        return result

    except Exception as e:
        print(f'fit_inhibitor_model ({model}): fitting failed – {e}')
        return None


# ── Fit table helper ──────────────────────────────────────────────────────────────

def _elongation_rate(fr, gene_length_aa, drug_diffusion_min):
    """Elongation rate in aa/sec from a fit_inhibitor_model result dict."""
    dwell = fr['t_runoff'] - drug_diffusion_min
    if dwell <= 0 or not np.isfinite(dwell):
        return np.nan
    return gene_length_aa / (dwell * 60.0)


def _print_fit_table(fit_results, fit_start_idx, fit_model, r2_threshold=0.95,
                     gene_length_effective=None, drug_diffusion_time_min=1.0):
    """Print a fit comparison table.

    Marks the last endpoint where R² ≥ r2_threshold as RECOMMENDED.
    This is the fit with the most data points that still meets the
    R² quality criterion.
    """
    valid = [r for r in fit_results if r['fit_result'] is not None]
    if not valid:
        print('  (no valid fits to tabulate)')
        return

    # Last endpoint with R² >= threshold → most data points still qualifying
    last_good_r2 = None
    for r in valid:
        if r['fit_result']['R2'] >= r2_threshold:
            last_good_r2 = r   # keep updating so the last qualifying row wins

    show_elong = gene_length_effective is not None
    w = 76 if show_elong else 60
    print(f'\n{"─"*w}')
    print(f' Fit Comparison  (start={fit_start_idx}, model={fit_model}, R²≥{r2_threshold})')
    print(f'{"─"*w}')
    hdr = f' {"#":>3}  {"end_idx":>7}  {"n_pts":>5}  {"k":>2}  {"χ²_red":>8}  {"R²":>7}'
    if show_elong:
        hdr += f'  {"ke(aa/s)":>12}'
    print(hdr)
    print(f'{"─"*w}')
    for i, r in enumerate(fit_results, 1):
        fr = r['fit_result']
        if fr is None:
            print(f' {i:>3}  {r["fit_end_idx"]:>7}  --- fit failed ---')
            continue
        n     = fr['n_data']
        k     = n - fr['dof']
        chi2r = fr['chi2_reduced']
        if last_good_r2 is not None and r['fit_end_idx'] == last_good_r2['fit_end_idx']:
            marker = f'  ← RECOMMENDED'
        else:
            marker = ''
        row = f' {i:>3}  {r["fit_end_idx"]:>7}  {n:>5}  {k:>2}  {chi2r:>8.4f}  {fr["R2"]:>7.3f}'
        if show_elong:
            er = r.get('elong_rate', np.nan)
            row += f'  {er:>12.3f}' if np.isfinite(er) else f'  {"---":>12}'
        row += marker
        print(row)
    print(f'{"─"*w}')
    if last_good_r2 is None:
        print(f' No endpoint found with R² ≥ {r2_threshold}. Lower r2_threshold.')
    print(f'{"─"*w}\n')


def plot_inhibitor(full_frames, intensities_normalized, inhibitor_frame_index,
                   results_folder=None, plot_name='HT', list_param=None,
                   responding_indices=None, figsize=(6, 3), time_array_min=None,
                   mean_intensity_ssa_inh=None, err_intensity_ssa_inh=None,
                   use_sem=True, show_individual_trajectories=True,
                   ylims=(0, 1.5), xlims=None,
                   y_label='Norm. Intensity',
                   treatment_label='Inhibitor', show_treatment_line=True,
                   # ── New fitting parameters ──
                   fit_model=None, fit_start_idx=None, fit_end_idx=None,
                   show_fit=True, show_runoff_time=True,
                   colors=None,
                   runoff_fraction=0.95,
                   remove_background_intensity=False,
                   background_frames=10,
                   show_background_line=False,
                   show_zero_y_axis=False,
                   fit_end_range=None, r2_threshold=0.95,
                   gene_length_effective=None, drug_diffusion_time_min=1.0,
                   frame_interval_sec=60):
    """Plot inhibitor run-off data with optional model fit.

    Args:
        full_frames: Time array (e.g., minutes, recentered so 0 = inhibitor).
        intensities_normalized: 2D array (n_cells × n_frames) of normalized
            intensities.
        inhibitor_frame_index: Frame index at which treatment starts.
        fit_model: Model to fit: 'linear', 'exponential', 'heaviside', or
            None (no fit).
        fit_start_idx: Start index for fitting range. Defaults to
            inhibitor_frame_index (t=0).
        fit_end_idx: End index for fitting range (inclusive). Defaults to
            last frame.
        show_fit: If True (default), overlay the fitted curve on the plot.
        show_runoff_time: If True (default), draw vertical lines for t½ and
            τ_runoff.
        runoff_fraction: Fraction of total decay for run-off time definition
            (default 0.95).
        colors: List of colors for the trajectories. If None, default colors
            are used.
        remove_background_intensity: If True, subtract the background intensity
            (estimated from the last `background_frames` frames within the
            xlims range) and rescale each trace to [0, 1] using the
            pre-treatment mean. Default False.
        background_frames: Number of frames at the end of the experiment (or
            xlims window) used to estimate background intensity. Default 10.
        show_background_line: If True, draw a horizontal dashed line at the
            estimated background intensity level. Default False.
        show_zero_y_axis: If True, draw a horizontal dashed line at y = 0.
        fit_end_range: When provided as (start, end), performs a sweep: fits
            the model from fit_start_idx to each endpoint in
            range(start, end+1). Only meaningful with fit_model='linear'.
            Endpoint values are automatically clamped to the available data
            length. Ignored when None (default).

    Returns:
        When fit_end_range is None: the fit result dict, or None.
        When fit_end_range is set: list of dicts per endpoint containing
        'fit_end_idx' and 'fit_result'.
    """
    if colors is None:
        colors = ['blue']
    if not isinstance(colors, list):
        colors = [colors]

    if results_folder is None:
        results_folder = current_dir / 'results_HT'
        results_folder.mkdir(exist_ok=True)

    if intensities_normalized is None or len(intensities_normalized) == 0:
        print('No data to plot.')
        return None

    if responding_indices is None:
        responding_indices = list(range(len(intensities_normalized)))
    if len(responding_indices) == 0:
        print('No responding cells to plot.')
        return None

    # ── Background removal and 0-1 rescaling ────────────────────────
    _bg_raw_value = None  # store for optional dashed line
    if remove_background_intensity:
        # Determine the end frame index from xlims or use all data
        if xlims is not None:
            end_mask = full_frames <= xlims[1]
            end_idx = int(np.sum(end_mask))
        else:
            end_idx = intensities_normalized.shape[1]
        bg_start = max(0, end_idx - background_frames)
        # Estimate background from the MEAN trajectory (robust to per-trace noise)
        resp_data = intensities_normalized[responding_indices, :]
        mean_for_bg = np.nanmean(resp_data[:, bg_start:end_idx])
        _bg_raw_value = mean_for_bg  # save original background value
        intensities_normalized = intensities_normalized - mean_for_bg
        # Rescale so the mean pre-treatment intensity = 1
        mean_pre = np.nanmean(resp_data[:, :inhibitor_frame_index] - mean_for_bg)
        if mean_pre != 0:
            intensities_normalized = intensities_normalized / mean_pre

    # Mean ± error (NaN-safe for artifact-removed frames) — needed by both paths
    if responding_indices:
        mean_trajectory = np.nanmean(intensities_normalized[responding_indices, :], axis=0)
        std_trajectory = np.nanstd(intensities_normalized[responding_indices, :], axis=0)
        if use_sem:
            n_valid = np.sum(np.isfinite(intensities_normalized[responding_indices, :]), axis=0)
            n_valid = np.maximum(n_valid, 1)
            err_trajectory = std_trajectory / np.sqrt(n_valid)
        else:
            err_trajectory = std_trajectory

    # ── Compute basal value for linear_extrapolated ───────────────────
    _basal_value = None
    if fit_model == 'linear_extrapolated':
        if remove_background_intensity:
            # After bg removal + rescaling, basal is 0
            _basal_value = 0.0
        else:
            # Use last background_frames of the normalized mean trajectory
            if xlims is not None:
                _end_mask = full_frames <= xlims[1]
                _basal_end_idx = int(np.sum(_end_mask))
            else:
                _basal_end_idx = len(mean_trajectory)
            _basal_bg_start = max(0, _basal_end_idx - background_frames)
            _basal_value = float(np.nanmean(mean_trajectory[_basal_bg_start:_basal_end_idx]))

    # ── AIC sweep (no main plot generated) ───────────────────────────────────
    if fit_end_range is not None:
        if fit_model not in ('linear', 'linear_extrapolated'):
            print(f'Warning: fit_end_range is only supported with fit_model="linear" '
                  f'or "linear_extrapolated". '
                  f'Ignoring fit_end_range (got fit_model={fit_model!r}).')
            return None

        start = fit_start_idx if fit_start_idx is not None else inhibitor_frame_index
        n_frames = len(full_frames)
        end_min = min(fit_end_range[0], n_frames - 1)
        end_max = min(fit_end_range[1], n_frames - 1)

        fit_results = []
        for end_idx in range(end_min, end_max + 1):
            if end_idx <= start:
                continue
            fr = fit_inhibitor_model(
                full_frames, mean_trajectory,
                err_data=err_trajectory,
                model=fit_model,
                fit_start_idx=start,
                fit_end_idx=end_idx,
                runoff_fraction=runoff_fraction,
                basal_value=_basal_value,
            )
            _er = (_elongation_rate(fr, gene_length_effective, drug_diffusion_time_min)
                   if (fr is not None and gene_length_effective is not None) else np.nan)
            fit_results.append({'fit_end_idx': end_idx, 'fit_result': fr, 'elong_rate': _er})

        # Per-endpoint plots
        for r in fit_results:
            end_idx = r['fit_end_idx']
            fr = r['fit_result']
            fig2, ax2 = plt.subplots(figsize=figsize, facecolor='white')
            ax2.set_facecolor('white')
            ax2.plot(full_frames, mean_trajectory, 'o-',
                     color=colors[0], linewidth=1, label='Experimental (mean)', markersize=6)
            ax2.fill_between(full_frames,
                             mean_trajectory - err_trajectory,
                             mean_trajectory + err_trajectory,
                             color=colors[0], alpha=0.2)
            if show_treatment_line:
                ax2.axvline(x=0, color='black', linestyle='--', linewidth=1,
                            label=f'{treatment_label} Treatment')
            if fr is not None:
                # Crop to the actual fitting window [start : end_idx+1]
                fit_x = full_frames[start:end_idx + 1]
                fit_y = fr['fitted_curve'][start:end_idx + 1]
                positive_mask = fit_y > 0
                if np.any(positive_mask):
                    last_pos = np.where(positive_mask)[0][-1] + 1
                    fit_x = fit_x[:last_pos]
                    fit_y = fit_y[:last_pos]
                else:
                    fit_x = fit_x[:0]
                    fit_y = fit_y[:0]
                if len(fit_x) > 0:
                    _fit_label = f'Lin.Extrap. Fit' if fit_model == 'linear_extrapolated' else f'Linear Fit'
                    ax2.plot(fit_x, fit_y, '-', color='red', linewidth=1.5,
                             label=f'{_fit_label} ({start}, {end_idx})\nR²={fr["R2"]:.3f}')

                # Dashed red extrapolation line for linear_extrapolated
                if fr['model'] == 'linear_extrapolated' and len(fit_x) > 0:
                    _t_ro = fr['t_runoff']
                    if np.isfinite(_t_ro) and _t_ro > fit_x[-1]:
                        _extrap_x = np.linspace(fit_x[-1], _t_ro, 50)
                        _extrap_y = _linear_model(_extrap_x, fr['params']['a (slope)'],
                                                  fr['params']['b (intercept)'])
                        ax2.plot(_extrap_x, _extrap_y, '--', color='red', linewidth=1.5,
                                 label='_nolegend_')

                if show_runoff_time and np.isfinite(fr['t_half']) and np.isfinite(fr['t_runoff']):
                    _er = r.get('elong_rate', np.nan)
                    _elong_lbl = f'  ({_er:.2f} aa/s)' if np.isfinite(_er) else ''

                    if fr['model'] != 'linear_extrapolated':
                        ax2.axvline(x=fr['t_half'], color='green', linestyle='--', linewidth=1,
                                    label=fr'$t_{{1/2}}$ ~ {fr["t_half"]:.1f} min')
                    ax2.axvline(x=fr['t_runoff'], color='orange', linestyle='--', linewidth=1,
                                label=fr'$\tau_{{runoff}}$ ~ {fr["t_runoff"]:.1f} min{_elong_lbl}')

            if show_background_line and _bg_raw_value is not None:
                _bg_y = 0 if remove_background_intensity else _bg_raw_value
                ax2.axhline(y=_bg_y, color='gray', linestyle=':', linewidth=1,
                            label='_nolegend_')

            ax2.set_xlabel('Time (min)', fontdict={'size': 16, 'color': 'black'})
            ax2.set_ylabel(y_label, fontdict={'size': 16, 'color': 'black'})
            ax2.tick_params(axis='both', which='major', labelsize=16,
                            labelcolor='black', colors='black')
            for spine in ax2.spines.values():
                spine.set_color('black')
                spine.set_linewidth(1.5)
            ax2.set_ylim(ylims)
            if xlims is not None:
                ax2.set_xlim(xlims)
            fig2.tight_layout()
            leg2 = ax2.legend(fontsize=10, loc='center left',
                              bbox_to_anchor=(1.02, 0.5),
                              framealpha=0.9, edgecolor='black')
            fname = f'HT_{plot_name}_end_{end_idx}'
            fig2.savefig(results_folder / (fname + '.png'), dpi=600,
                         bbox_extra_artists=(leg2,), bbox_inches='tight')
            fig2.savefig(results_folder / (fname + '.svg'), dpi=600,
                         bbox_extra_artists=(leg2,), bbox_inches='tight')
            plt.show()

        _print_fit_table(fit_results, start, fit_model, r2_threshold=r2_threshold,
                         gene_length_effective=gene_length_effective,
                         drug_diffusion_time_min=drug_diffusion_time_min)
        return fit_results

    # ── Main plot (only when fit_end_range is None) ───────────────────────────
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    ax.set_facecolor('white')

    # Individual trajectories
    if show_individual_trajectories:
        if responding_indices:
            for i in responding_indices:
                ax.plot(full_frames, intensities_normalized[i],
                        linestyle='-', color='dimgray', linewidth=0.2)

    ax.plot(full_frames, mean_trajectory, 'o-',
            color=colors[0], linewidth=1, label='Experimental (mean)', markersize=6)
    ax.fill_between(full_frames,
                    mean_trajectory - err_trajectory,
                    mean_trajectory + err_trajectory,
                    color=colors[0], alpha=0.2)

    # TASEP simulation overlay (if provided)
    if mean_intensity_ssa_inh is not None and err_intensity_ssa_inh is not None:
        _has_params = (list_param is not None
                       and list_param[0] is not None
                       and list_param[1] is not None)
        legend_label_sim = (fr'Model Fit ($k_e$={np.round(list_param[1],1)}, $k_i$={np.round(list_param[0],3)})'
                        if _has_params else 'Simulation')
        _treatment_time_min = inhibitor_frame_index * frame_interval_sec / 60.0
        ax.plot(time_array_min - _treatment_time_min, mean_intensity_ssa_inh, '-', color='red', linewidth=3, label=legend_label_sim)
        ax.fill_between(time_array_min - _treatment_time_min, mean_intensity_ssa_inh - err_intensity_ssa_inh,
                        mean_intensity_ssa_inh + err_intensity_ssa_inh, color='red', alpha=0.1)

    # ── Model fit ────────────────────────────────────────────────────────
    fit_result = None
    if fit_model is not None and fit_end_range is None:
        # Default fit range: from inhibitor application to end
        start = fit_start_idx if fit_start_idx is not None else inhibitor_frame_index
        end = fit_end_idx  # None → full array end (handled inside fit_inhibitor_model)

        fit_result = fit_inhibitor_model(
            full_frames, mean_trajectory,
            err_data=err_trajectory,
            model=fit_model,
            fit_start_idx=start,
            fit_end_idx=end,
            runoff_fraction=runoff_fraction,
            basal_value=_basal_value,
        )

        if fit_result is not None:
            _end_str = fit_end_idx if fit_end_idx is not None else "end"
            model_labels = {'linear': f'Linear Fit ({start}, {_end_str})', 'exponential': 'Exponential Fit',
                            'heaviside': 'Heaviside Fit',
                            'linear_extrapolated': f'Lin.Extrap. Fit ({start}, {_end_str})'}
            label = model_labels.get(fit_result['model'], 'Fit')

            t_half = fit_result['t_half']
            t_runoff = fit_result['t_runoff']

            if show_fit:
                _end = (fit_end_idx + 1) if fit_end_idx is not None else len(full_frames)
                fit_x = full_frames[start:_end]
                fit_y = fit_result['fitted_curve'][start:_end]
                if fit_result['model'] in ('linear', 'linear_extrapolated'):
                    positive_mask = fit_y > 0
                    if np.any(positive_mask):
                        last_pos = np.where(positive_mask)[0][-1] + 1
                        fit_x = fit_x[:last_pos]
                        fit_y = fit_y[:last_pos]
                    else:
                        fit_x = fit_x[:0]
                        fit_y = fit_y[:0]
                if len(fit_x) > 0:
                    ax.plot(fit_x, fit_y, '-', color='red', linewidth=1.5, label=label + f'\nR²={fit_result["R2"]:.3f}')

                # Dashed red extrapolation line for linear_extrapolated
                if fit_result['model'] == 'linear_extrapolated' and len(fit_x) > 0:
                    if np.isfinite(t_runoff) and t_runoff > fit_x[-1]:
                        _extrap_x = np.linspace(fit_x[-1], t_runoff, 50)
                        _extrap_y = _linear_model(_extrap_x, fit_result['params']['a (slope)'],
                                                  fit_result['params']['b (intercept)'])
                        ax.plot(_extrap_x, _extrap_y, '--', color='red', linewidth=1.5,
                                label='_nolegend_')

            _er = (_elongation_rate(fit_result, gene_length_effective, drug_diffusion_time_min)
                   if gene_length_effective is not None else np.nan)

            if show_runoff_time:
                if fit_result['model'] != 'linear_extrapolated':
                    ax.axvline(x=t_half, color='green', linestyle='--', linewidth=1,
                               label=fr'$t_{{1/2}}$ ~ {t_half:.1f} min')
                _elong_lbl = f'  ({_er:.2f} aa/s)' if np.isfinite(_er) else ''
                ax.axvline(x=t_runoff, color='orange', linestyle='--', linewidth=1,
                           label=fr'$\tau_{{ro}}$ {t_runoff:.1f} min{_elong_lbl}')

            # Print fitted parameters
            print(f'── {label} ──')
            for k, v in fit_result['params'].items():
                print(f'  {k}: {v:.4f}')
            print(f'  t½:      {fit_result["t_half"]:.2f} min')
            if fit_result['model'] == 'linear_extrapolated':
                print(f'  τ_runoff (basal crossing): {fit_result["t_runoff"]:.2f} min')
            else:
                print(f'  τ_runoff (2×t½): {fit_result["t_runoff"]:.2f} min')
            if gene_length_effective is not None and np.isfinite(_er):
                print(f'  ke:      {_er:.4f} aa/s')
            chi2r = fit_result['chi2_reduced']
            print(f'  χ²_red:  {chi2r:.4f}  (dof={fit_result["dof"]})')
            print(f'  R²:      {fit_result["R2"]:.3f}  (n={fit_result["n_data"]})')

    # Treatment line at t = 0
    if show_treatment_line:
        ax.axvline(x=0, color='black', linestyle='--', linewidth=1,
                    label=f'{treatment_label} Treatment')

    ax.set_xlabel("Time (min)", fontdict={'size': 16, 'color': 'black'})
    ax.set_ylabel(y_label, fontdict={'size': 16, 'color': 'black'})
    ax.tick_params(axis='both', which='major', labelsize=16, labelcolor='black', colors='black')

    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)

    # ── Optional dashed reference lines ───────────────────────────────
    if show_zero_y_axis:
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8,
                    label='_nolegend_')
    if show_background_line and _bg_raw_value is not None:
        _bg_y = 0 if remove_background_intensity else _bg_raw_value
        ax.axhline(y=_bg_y, color='gray', linestyle=':', linewidth=1,
                    label='_nolegend_')
    elif show_background_line and _bg_raw_value is None:
        if xlims is not None:
            end_mask = full_frames <= xlims[1]
            end_idx = int(np.sum(end_mask))
        else:
            end_idx = intensities_normalized.shape[1]
        bg_start = max(0, end_idx - background_frames)
        resp_data = intensities_normalized[responding_indices, :]
        _bg_display = np.nanmean(resp_data[:, bg_start:end_idx])
        ax.axhline(y=_bg_display, color='gray', linestyle=':', linewidth=1,
                    label='_nolegend_')

    ax.set_ylim(ylims)
    if xlims is not None:
        ax.set_xlim(xlims)
    fig.tight_layout()
    legend = ax.legend(fontsize=10, loc='center left', bbox_to_anchor=(1.02, 0.5),
                       framealpha=0.9, edgecolor='black')
    fig.savefig(results_folder / f'HT_{plot_name}.png', dpi=600,
                bbox_extra_artists=(legend,), bbox_inches='tight')
    fig.savefig(results_folder / f'HT_{plot_name}.svg', dpi=600,
                bbox_extra_artists=(legend,), bbox_inches='tight')

    plt.show()

    return fit_result


def plot_multiple_inhibitors(full_frames_list,
                                intensities_normalized_list,
                                inhibitor_frame_index,
                                results_folder=None,
                                plot_name='HT_multi',
                                responding_indices_list=None,
                                source_data_collector=None,
                                figsize=(6, 3),
                                colors=None,
                                legend_labels=None,
                                use_sem=True,
                                show_individual_trajectories=True,
                                ylims=(0, 1.5),
                                xlims=None,
                                time_unit='min',
                                y_label='Norm. Intensity',
                                treatment_label='Inhibitor',
                                show_treatment_line=True,
                                show_n_cells=False,
                                n_cells_override=None,
                                # ── Fitting parameters ──
                                fit_model=None,
                                fit_start_idx=None,
                                fit_end_idx=None,
                                show_fit=True,
                                show_runoff_time=True,
                                show_half_life_lines=False,
                                show_half_life_table=False,
                                half_life_target=None,
                                half_life_exclude_labels=(),
                                runoff_fraction=0.95,
                                remove_background_intensity=False,
                                background_frames=10,
                                show_background_line=False,
                                show_zero_y_axis=False,
                                fit_end_range=None, r2_threshold=0.95,
                                gene_length_effective=None, drug_diffusion_time_min=1.0,
                                filename_prefix='HT_',
                                error_style='band',
                                capsize=4,
                                verbose=True,
                                connect_points=False,
                                font_family='Arial',
                                axis_label_size=20,
                                tick_label_size=18,
                                legend_font_size=14,
                                marker_size=8,
                                data_line_width=2.5,
                                error_bar_line_width=2.0,
                                error_bar_cap_thickness=1.5,
                                fit_line_width=2.5,
                                legend_position='right'):
    """Plot multiple inhibitor datasets on the same axes with optional model fits.

    Args:
        full_frames_list: Either a single 1D array of frame times (applies to
            all datasets) or a list of 1D arrays, one per dataset.
        intensities_normalized_list: List of (n_cells × n_frames) arrays of
            normalized intensities.
        inhibitor_frame_index: Frame index at which inhibitor treatment starts.
        results_folder: Where to save the figure (will be created if needed).
        plot_name: Filename suffix for the saved figure.
        filename_prefix: Prefix added to saved figure filenames. Defaults to
            ``HT_`` for backward compatibility.
        responding_indices_list: Per-dataset lists of cell indices to include.
            Defaults to all.
        source_data_collector: Optional list that receives one dictionary per
            plotted dataset containing the exact time, mean, error, sample-size,
            and dense fitted-curve arrays used to draw the figure. This is for
            publication source-data export and does not alter calculations.
        figsize: Figure size tuple.
        colors: Matplotlib color codes for each dataset.
        legend_labels: Text labels for each dataset's mean trace.
        show_n_cells: If True, append ``(n=<count>)`` to each legend label
            showing the number of cells (trajectories) in that dataset.
        n_cells_override: Optional list of integers, one per dataset. When
            provided, these counts are used instead of the row count of the
            data matrix.  Useful in ``comparison='average'`` mode where rows
            are repetition means, not individual cells.
        use_sem: If True, error bands show SEM; else SD.
        show_individual_trajectories: If True, plot individual cell traces.
        ylims: (ymin, ymax) for the plot.
        xlims: (xmin, xmax) for the plot, expressed in ``time_unit``. If
            None, the limits are auto-scaled.
        time_unit: Display unit for the x-axis: ``'min'`` (default) or
            ``'h'``/``'hours'``. Fitting and internal timing remain in
            minutes; when set to hours, pass ``xlims=(0, 4)`` for a
            four-hour display.
        fit_model: Model to fit per dataset: 'linear', 'exponential',
            'heaviside', or None.
        fit_start_idx: Start index for fitting range. Defaults to
            inhibitor_frame_index.
        fit_end_idx: End index for fitting range (inclusive). Defaults to
            last frame.
        show_fit: If True (default), overlay the fitted curve on the plot.
        show_runoff_time: If True (default), draw vertical lines for t½ and
            τ_runoff.
        show_half_life_lines: If True, draw condition-colored dashed guides
            for the fitted half-level and half-life. Intended for exponential
            fits.
        show_half_life_table: If True, print a table containing the fitted
            half-life and half-level for each dataset.
        half_life_target: Optional absolute normalized signal target. For
            example, ``0.5`` reports the time to reach 0.5 rather than the
            time to halfway between the fitted signal and plateau.
        half_life_exclude_labels: Dataset labels to omit from the half-life
            guides and table, such as ``('Control',)``.
        runoff_fraction: Fraction of total decay for run-off time definition
            (default 0.95).
        remove_background_intensity: If True, subtract background intensity
            and rescale each trace to [0, 1]. Default False.
        background_frames: Number of frames at the end used to estimate
            background intensity. Default 10.
        show_background_line: If True, draw a horizontal dashed line at the
            estimated background intensity level per dataset. Default False.
        show_zero_y_axis: If True, draw a horizontal dashed line at y = 0.
        error_style: ``'band'`` (default) renders mean ± error as a shaded
            ``fill_between`` region; ``'bar'`` renders discrete vertical
            error bars with caps at each time point via ``ax.errorbar``.
        capsize: Cap width in points for error bars when
            ``error_style='bar'``. Default 4.
        verbose: If True, print fitted parameters and goodness-of-fit
            diagnostics. Default True.
        connect_points: If True, connect the mean data points with straight
            line segments. This is independent of ``fit_model``. Default
            False.
        font_family: Font family used for axis labels, ticks, and legend.
        axis_label_size: Font size for the x- and y-axis labels.
        tick_label_size: Font size for x- and y-axis tick labels.
        legend_font_size: Font size for legend text.
        marker_size: Size of the plotted mean markers.
        data_line_width: Width of lines connecting the plotted mean markers.
        error_bar_line_width: Width of the error-bar lines.
        error_bar_cap_thickness: Thickness of the error-bar caps.
        fit_line_width: Width of fitted model curves.
        legend_position: Position of the legend: ``'right'`` (default, outside right),
            ``'top'`` (above axes), ``'inside_upper_right'`` / ``'upper right'``,
            ``'inside_lower_left'`` / ``'lower left'``, or any valid matplotlib ``loc`` string.
        fit_end_range: When provided as (start, end), performs a sweep per
            dataset. Only meaningful with fit_model='linear'. Endpoint values
            are automatically clamped to the available data length.

    Returns:
        When fit_end_range is None: list of fit result dicts per dataset.
        When fit_end_range is set: list of fit results lists per dataset,
        each element containing 'fit_end_idx' and 'fit_result'.
    """
    # Prepare output folder
    if results_folder is None:
        results_folder = current_dir / 'results_HT'
    results_folder.mkdir(parents=True, exist_ok=True)

    time_unit = str(time_unit).strip().lower()
    if time_unit in {'min', 'mins', 'minute', 'minutes'}:
        time_scale = 1.0
        time_axis_label = 'Time (min)'
    elif time_unit in {'h', 'hr', 'hrs', 'hour', 'hours'}:
        time_scale = 60.0
        time_axis_label = 'Time (h)'
    else:
        raise ValueError("time_unit must be 'min' or 'h'.")

    # Internal data and fit calculations remain in minutes.  xlims is in the
    # selected display unit so callers can use xlims=(0, 4) for hours.
    xlims_display = tuple(xlims) if xlims is not None else None
    xlims_min = (
        tuple(float(value) * time_scale for value in xlims_display)
        if xlims_display is not None else None
    )

    # Set up figure
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    ax.set_facecolor('white')

    # Default color cycle
    if colors is None:
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    fit_results = []
    half_life_rows = []
    half_life_exclude_labels = set(half_life_exclude_labels)
    if half_life_target is not None:
        half_life_target = float(half_life_target)
        if not np.isfinite(half_life_target):
            raise ValueError('half_life_target must be finite or None.')

    # Loop over each dataset
    for idx, intensities in enumerate(intensities_normalized_list):
        # Select frames array
        frames = (full_frames_list[idx]
                  if isinstance(full_frames_list, (list, tuple))
                  else full_frames_list)
        color = colors[idx % len(colors)]

        # Select responding cell indices
        resp_idx = (responding_indices_list[idx]
                    if (responding_indices_list and idx < len(responding_indices_list))
                    else list(range(intensities.shape[0])))

        # ── Background removal and 0-1 rescaling ────────────────
        _bg_raw_value = None
        frames_display = frames / time_scale
        if xlims_min is not None:
            end_mask = frames <= xlims_min[1]
            _end_idx = int(np.sum(end_mask))
        else:
            _end_idx = intensities.shape[1]
        _bg_start = max(0, _end_idx - background_frames)

        if remove_background_intensity:
            intensities = intensities.copy()
            # Estimate background from the MEAN trajectory (robust to per-trace noise)
            resp_data = intensities[resp_idx, :]
            mean_for_bg = np.nanmean(resp_data[:, _bg_start:_end_idx])
            _bg_raw_value = mean_for_bg
            intensities = intensities - mean_for_bg
            # Rescale so the mean pre-treatment intensity = 1
            mean_pre = np.nanmean(resp_data[:, :inhibitor_frame_index] - mean_for_bg)
            if mean_pre != 0:
                intensities = intensities / mean_pre

        if show_background_line:
            if _bg_raw_value is not None:
                # After rescaling, background maps to 0
                _lbl = legend_labels[idx] if legend_labels else f'Dataset {idx}'
                ax.axhline(y=0, color=color, linestyle=':', linewidth=1,
                            label=f'BG ({_lbl}: {_bg_raw_value:.1f} raw)')
            else:
                # No rescaling — compute and show at its original level
                resp_data = intensities[resp_idx, :]
                _bg_display = np.nanmean(resp_data[:, _bg_start:_end_idx])
                _lbl = legend_labels[idx] if legend_labels else f'Dataset {idx}'
                ax.axhline(y=_bg_display, color=color, linestyle=':', linewidth=1,
                            label=f'BG ({_lbl}) = {_bg_display:.1f}')

        # Plot individual trajectories
        if show_individual_trajectories:
            for i in resp_idx:
                ax.plot(frames_display, intensities[i],
                        linestyle='-', color=color,
                        linewidth=0.3, alpha=0.4)

        # Compute mean & error (NaN-safe for artifact-removed frames)
        data = intensities[resp_idx, :]
        mean_traj = np.nanmean(data, axis=0)
        if use_sem:
            n_valid = np.sum(np.isfinite(data), axis=0)
            squared_deviations = np.nansum(
                (data - mean_traj) ** 2,
                axis=0,
            )
            sample_variance = np.divide(
                squared_deviations,
                n_valid - 1,
                out=np.full_like(mean_traj, np.nan, dtype=float),
                where=n_valid > 1,
            )
            std_traj = np.sqrt(sample_variance)
            err_traj = np.divide(
                std_traj,
                np.sqrt(n_valid),
                out=np.full_like(std_traj, np.nan, dtype=float),
                where=n_valid > 1,
            )
        else:
            std_traj = np.nanstd(data, axis=0)
            err_traj = std_traj

        fit_x_export = np.asarray([], dtype=float)
        fit_y_export = np.asarray([], dtype=float)

        # Determine legend text
        condition_label_text = (legend_labels[idx]
                                if (legend_labels and idx < len(legend_labels))
                                else f'Dataset {idx+1}')
        label_text = condition_label_text
        if show_n_cells:
            _n_display = (n_cells_override[idx]
                          if (n_cells_override and idx < len(n_cells_override))
                          else len(resp_idx))
            label_text = f'{label_text} (n={_n_display})'

        # Plot mean ± error
        if error_style == 'bar':
            ax.errorbar(frames_display, mean_traj, yerr=err_traj,
                        # Plot observations as points only.  The smooth model
                        # curve below carries the temporal connection.
                        fmt='o-' if connect_points else 'o',
                        color=color, ecolor=color,
                        linestyle='-' if connect_points else 'None',
                        linewidth=data_line_width, elinewidth=error_bar_line_width,
                        markersize=marker_size,
                        capsize=capsize, capthick=error_bar_cap_thickness,
                        zorder=3,
                        label=label_text)
        else:
            ax.plot(frames_display, mean_traj,
                    'o-' if connect_points else 'o', color=color,
                    linewidth=data_line_width, markersize=marker_size, zorder=3,
                    label=label_text)
            ax.fill_between(frames_display,
                            mean_traj - err_traj,
                            mean_traj + err_traj,
                            color=color, alpha=0.2)

        # ── Compute basal value for linear_extrapolated ────────────
        _ds_basal_value = None
        if fit_model == 'linear_extrapolated':
            if remove_background_intensity:
                _ds_basal_value = 0.0
            else:
                if xlims_min is not None:
                    _basal_end = int(np.sum(frames <= xlims_min[1]))
                else:
                    _basal_end = len(mean_traj)
                _basal_start = max(0, _basal_end - background_frames)
                _ds_basal_value = float(np.nanmean(mean_traj[_basal_start:_basal_end]))

        # ── Model fit ────────────────────────────────────────────────
        if fit_model is not None and fit_end_range is None:
            start = fit_start_idx if fit_start_idx is not None else inhibitor_frame_index
            end = fit_end_idx

            fit_result = fit_inhibitor_model(
                frames, mean_traj,
                err_data=err_traj,
                model=fit_model,
                fit_start_idx=start,
                fit_end_idx=end,
                runoff_fraction=runoff_fraction,
                basal_value=_ds_basal_value,
            )

            if fit_result is not None:
                _end_str = fit_end_idx if fit_end_idx is not None else "end"
                model_labels = {'linear': f'Linear Fit ({start}, {_end_str})', 'exponential': 'Exponential Fit',
                                'heaviside': 'Heaviside Fit',
                                'linear_extrapolated': f'Lin.Extrap. Fit ({start}, {_end_str})'}
                fit_label = model_labels.get(fit_result['model'], 'Fit')
                t_half = fit_result['t_half']
                t_runoff = fit_result['t_runoff']

                # Half-life is measured from the first fitted post-treatment
                # time point. Use the observed mean at that point as S(t0),
                # then calculate C + [S(t0) - C] / 2.
                if (fit_result['model'] == 'exponential'
                        and condition_label_text not in half_life_exclude_labels):
                    _fit_params = fit_result['params']
                    _t0_min = float(frames[start])
                    _s0 = float(mean_traj[start])
                    if not np.isfinite(_s0):
                        _s0 = float(_exponential_model(
                            _t0_min,
                            _fit_params['A (amplitude)'],
                            _fit_params['tau (time constant)'],
                            _fit_params['C (baseline)'],
                        ))
                    _baseline = float(_fit_params['C (baseline)'])
                    if half_life_target is None:
                        _half_level = float(_baseline + 0.5 * (_s0 - _baseline))
                        _half_duration = float(t_half)
                    else:
                        _half_level = float(half_life_target)
                        _target_component = _half_level - _baseline
                        _initial_component = _s0 - _baseline
                        if (_initial_component > 0
                                and _target_component > 0
                                and _target_component < _initial_component):
                            _half_duration = float(
                                _fit_params['tau (time constant)']
                                * np.log(_initial_component / _target_component)
                            )
                        else:
                            _half_duration = np.nan
                    _half_time_min = (_t0_min + _half_duration
                                      if np.isfinite(_half_duration) else np.nan)
                    half_life_rows.append({
                        'Condition': condition_label_text,
                        'Half time (min)': _half_duration,
                        'Half time (h)': (_half_duration / 60.0
                                         if np.isfinite(_half_duration)
                                         else np.nan),
                    })

                    if show_half_life_lines and np.isfinite(_half_time_min):
                        _fit_end_index = min(
                            fit_end_idx if fit_end_idx is not None else len(frames) - 1,
                            len(frames) - 1,
                        )
                        _line_end = frames[_fit_end_index]
                        ax.plot(
                            [_t0_min / time_scale, _line_end / time_scale],
                            [_half_level, _half_level],
                            linestyle='--', color=color, linewidth=1,
                            alpha=0.65, label='_nolegend_', zorder=0,
                        )
                        ax.vlines(
                            x=_half_time_min / time_scale,
                            ymin=0,
                            ymax=_half_level,
                            linestyle='--', color=color, linewidth=1,
                            alpha=0.65, label='_nolegend_', zorder=0,
                        )

                if show_fit:
                    _end = (fit_end_idx + 1) if fit_end_idx is not None else len(frames)
                    fit_x_data = frames[start:_end]
                    if len(fit_x_data) == 0:
                        fit_x = fit_x_data
                        fit_y = fit_x_data
                    else:
                        # Evaluate the fitted model on a dense grid so the
                        # curve looks continuous rather than like another
                        # line connecting the measured time points.
                        fit_x = np.linspace(fit_x_data[0], fit_x_data[-1], 300)
                        fit_params = fit_result['params']
                        if fit_result['model'] == 'exponential':
                            fit_y = _exponential_model(
                                fit_x,
                                fit_params['A (amplitude)'],
                                fit_params['tau (time constant)'],
                                fit_params['C (baseline)'],
                            )
                        elif fit_result['model'] == 'heaviside':
                            fit_y = _heaviside_model(
                                fit_x,
                                fit_params['A (amplitude)'],
                                fit_params['T (dwell/run-off time)'],
                                fit_params['C (baseline)'],
                            )
                        elif fit_result['model'] == 'linear':
                            fit_y = _linear_model(
                                fit_x,
                                fit_params['a (slope)'],
                                fit_params['b (intercept)'],
                            )
                        else:
                            # Keep the fallback for any future model types.
                            fit_y = np.interp(
                                fit_x, fit_x_data,
                                fit_result['fitted_curve'][start:_end],
                            )

                    if fit_result['model'] in ('linear', 'linear_extrapolated'):
                        positive_mask = fit_y > 0
                        if np.any(positive_mask):
                            last_pos = np.where(positive_mask)[0][-1] + 1
                            fit_x = fit_x[:last_pos]
                            fit_y = fit_y[:last_pos]
                        else:
                            fit_x = fit_x[:0]
                            fit_y = fit_y[:0]
                    if len(fit_x) > 0:
                        ax.plot(fit_x / time_scale, fit_y, '-',
                                color=color, linewidth=fit_line_width, alpha=0.9,
                                solid_capstyle='round', zorder=1,
                                label='_nolegend_')
                        fit_x_export = np.asarray(fit_x / time_scale, dtype=float).copy()
                        fit_y_export = np.asarray(fit_y, dtype=float).copy()

                    # Same-color extrapolation for linear_extrapolated
                    if fit_result['model'] == 'linear_extrapolated' and len(fit_x) > 0:
                        if np.isfinite(t_runoff) and t_runoff > fit_x[-1]:
                            _extrap_x = np.linspace(fit_x[-1], t_runoff, 50)
                            _extrap_y = _linear_model(_extrap_x, fit_result['params']['a (slope)'],
                                                      fit_result['params']['b (intercept)'])
                            ax.plot(_extrap_x / time_scale, _extrap_y, '--', color=color,
                                    linewidth=fit_line_width,
                                    label='_nolegend_')

                if show_runoff_time:
                    _ro_label = 'basal' if fit_result['model'] == 'linear_extrapolated' else '2×t½'
                    if fit_result['model'] != 'linear_extrapolated':
                        ax.axvline(x=t_half / time_scale, color=color, linestyle=':', linewidth=1,
                                   label=f'{label_text} t½ ~ {t_half / time_scale:.1f} {"h" if time_scale == 60 else "min"}')
                    ax.axvline(x=t_runoff / time_scale, color=color, linestyle='--', linewidth=1,
                               label=f'{label_text} τ ({_ro_label}) ~ {t_runoff / time_scale:.1f} {"h" if time_scale == 60 else "min"}')

                if verbose:
                    # Print fitted parameters
                    print(f'── {label_text}: {fit_label} ──')
                    for k, v in fit_result['params'].items():
                        print(f'  {k}: {v:.4f}')
                    print(f'  t½:      {fit_result["t_half"]:.2f} min')
                    if fit_result['model'] == 'linear_extrapolated':
                        print(f'  τ_runoff (basal crossing): {fit_result["t_runoff"]:.2f} min')
                    else:
                        print(f'  τ_runoff (2×t½): {fit_result["t_runoff"]:.2f} min')
                    chi2r = fit_result['chi2_reduced']
                    print(f'  χ²_red:  {chi2r:.4f}  (dof={fit_result["dof"]})')
                    print(f'  R²:      {fit_result["R2"]:.3f}  (n={fit_result["n_data"]})')

            fit_results.append(fit_result)
        else:
            fit_results.append(None)

        if source_data_collector is not None:
            n_cells_total = (
                n_cells_override[idx]
                if n_cells_override and idx < len(n_cells_override)
                else len(resp_idx)
            )
            source_data_collector.append({
                'dataset_index': idx,
                'condition_label': condition_label_text,
                'time': np.asarray(frames_display, dtype=float).copy(),
                'mean': np.asarray(mean_traj, dtype=float).copy(),
                'error': np.asarray(err_traj, dtype=float).copy(),
                'error_type': 'sem' if use_sem else 'sd',
                'n_cells_total': int(n_cells_total),
                'n_independent_units': int(len(resp_idx)),
                'fit_time': fit_x_export,
                'fit_value': fit_y_export,
            })

    # Plot treatment line at zero
    if show_treatment_line:
        ax.axvline(x=0, color='black', linestyle='--', linewidth=1,
                    label=treatment_label)

    # Optional horizontal reference lines
    if show_zero_y_axis:
        ax.axhline(y=0, color='black', linestyle='--', linewidth=0.8,
                    label='y = 0')

    # Styling
    ax.set_xlabel(time_axis_label, fontsize=axis_label_size,
                  fontname=font_family, color='black')
    ax.set_ylabel(y_label, fontsize=axis_label_size,
                  fontname=font_family, color='black')
    ax.tick_params(axis='both', which='major', labelsize=tick_label_size,
                   labelcolor='black', colors='black')
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontname(font_family)
        tick_label.set_color('black')
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)
    ax.set_ylim(*ylims)
    if xlims_display is not None:
        ax.set_xlim(*xlims_display)

    if legend_position == 'top':
        legend = ax.legend(
            fontsize=legend_font_size,
            loc='lower center',
            bbox_to_anchor=(0.5, 1.02),
            ncol=len(legend_labels),
            frameon=False,
            handlelength=2.5,
            columnspacing=1.5,
        )
    elif legend_position in {'inside_lower_left', 'lower left'}:
        legend = ax.legend(
            fontsize=legend_font_size,
            loc='lower left',
            frameon=False,
        )
    elif legend_position in {'inside_upper_right', 'inside_top_right', 'upper right', 'top right'}:
        legend = ax.legend(
            fontsize=legend_font_size,
            loc='upper right',
            frameon=False,
        )
    elif legend_position in {'right', 'outside_right'}:
        legend = ax.legend(fontsize=legend_font_size, loc='center left',
                           bbox_to_anchor=(1.02, 0.5), frameon=False)
    else:
        legend = ax.legend(
            fontsize=legend_font_size,
            loc=legend_position,
            frameon=False,
        )
    for legend_text in legend.get_texts():
        legend_text.set_fontname(font_family)
        legend_text.set_color('black')
    fig.tight_layout()

    if show_half_life_table and half_life_rows:
        _half_life_table = pd.DataFrame(half_life_rows)
        print('\nHalf-life summary')
        print(_half_life_table.to_string(
            index=False,
            formatters={
                'Half time (min)': '{:.2f}'.format,
                'Half time (h)': '{:.2f}'.format,
            },
        ))

    # Save & show
    fig.savefig(results_folder / f'{filename_prefix}{plot_name}.png', dpi=600,
                bbox_extra_artists=(legend,), bbox_inches='tight')
    fig.savefig(results_folder / f'{filename_prefix}{plot_name}.svg', dpi=600,
                bbox_extra_artists=(legend,), bbox_inches='tight')
    plt.show()

    # ── AIC sweep (per dataset) ───────────────────────────────────────────────
    if fit_end_range is not None:
        if fit_model not in ('linear', 'linear_extrapolated'):
            print(f'Warning: fit_end_range is only supported with fit_model="linear" '
                  f'or "linear_extrapolated". '
                  f'Ignoring fit_end_range (got fit_model={fit_model!r}).')
            return fit_results

        start = fit_start_idx if fit_start_idx is not None else inhibitor_frame_index
        all_fit_results = []

        for ds_idx, intensities in enumerate(intensities_normalized_list):
            frames = (full_frames_list[ds_idx]
                      if isinstance(full_frames_list, (list, tuple))
                      else full_frames_list)
            color = colors[ds_idx % len(colors)]
            resp_idx = (responding_indices_list[ds_idx]
                        if (responding_indices_list
                            and ds_idx < len(responding_indices_list))
                        else list(range(intensities.shape[0])))
            ds_label = (legend_labels[ds_idx]
                        if (legend_labels and ds_idx < len(legend_labels))
                        else f'Dataset {ds_idx+1}')

            # Re-compute mean/err for this dataset (mirrors the main loop above)
            _int = intensities.copy()
            _bg_raw_value = None  # per-dataset; avoids leakage from the main loop
            if remove_background_intensity:
                if xlims_min is not None:
                    _end_idx = int(np.sum(frames <= xlims_min[1]))
                else:
                    _end_idx = _int.shape[1]
                _bg_start = max(0, _end_idx - background_frames)
                _resp_data = _int[resp_idx, :]
                _mean_bg = np.nanmean(_resp_data[:, _bg_start:_end_idx])
                _bg_raw_value = _mean_bg
                _int = _int - _mean_bg
                _mean_pre = np.nanmean(_resp_data[:, :inhibitor_frame_index] - _mean_bg)
                if _mean_pre != 0:
                    _int = _int / _mean_pre

            _data = _int[resp_idx, :]
            _mean_traj = np.nanmean(_data, axis=0)
            _std_traj = np.nanstd(_data, axis=0)
            if use_sem:
                _n_valid = np.maximum(np.sum(np.isfinite(_data), axis=0), 1)
                _err_traj = _std_traj / np.sqrt(_n_valid)
            else:
                _err_traj = _std_traj

            # Compute basal value for linear_extrapolated (per dataset)
            _sweep_basal = None
            if fit_model == 'linear_extrapolated':
                if remove_background_intensity:
                    _sweep_basal = 0.0
                else:
                    if xlims_min is not None:
                        _sb_end = int(np.sum(frames <= xlims_min[1]))
                    else:
                        _sb_end = len(_mean_traj)
                    _sb_start = max(0, _sb_end - background_frames)
                    _sweep_basal = float(np.nanmean(_mean_traj[_sb_start:_sb_end]))

            n_frames = len(frames)
            end_min = min(fit_end_range[0], n_frames - 1)
            end_max = min(fit_end_range[1], n_frames - 1)

            fit_results = []
            for end_idx in range(end_min, end_max + 1):
                if end_idx <= start:
                    continue
                fr = fit_inhibitor_model(
                    frames, _mean_traj,
                    err_data=_err_traj,
                    model=fit_model,
                    fit_start_idx=start,
                    fit_end_idx=end_idx,
                    runoff_fraction=runoff_fraction,
                    basal_value=_sweep_basal,
                )
                _er = (_elongation_rate(fr, gene_length_effective, drug_diffusion_time_min)
                       if (fr is not None and gene_length_effective is not None) else np.nan)
                fit_results.append({'fit_end_idx': end_idx, 'fit_result': fr, 'elong_rate': _er})

            # Per-endpoint plots for this dataset
            for r in fit_results:
                end_idx = r['fit_end_idx']
                fr = r['fit_result']
                fig2, ax2 = plt.subplots(figsize=figsize, facecolor='white')
                ax2.set_facecolor('white')
                frames_display = frames / time_scale
                ax2.plot(frames_display, _mean_traj, 'o-', color=color,
                         linewidth=1, label=ds_label, markersize=6)
                ax2.fill_between(frames_display,
                                 _mean_traj - _err_traj,
                                 _mean_traj + _err_traj,
                                 color=color, alpha=0.2)
                if show_treatment_line:
                    ax2.axvline(x=0, color='black', linestyle='--', linewidth=1,
                                label=f'{treatment_label} Treatment')
                if fr is not None:
                    # Crop to the actual fitting window [start : end_idx+1]
                    fit_x = frames[start:end_idx + 1]
                    fit_y = fr['fitted_curve'][start:end_idx + 1]
                    positive_mask = fit_y > 0
                    if np.any(positive_mask):
                        lp = np.where(positive_mask)[0][-1] + 1
                        fit_x = fit_x[:lp]
                        fit_y = fit_y[:lp]
                    else:
                        fit_x = fit_x[:0]
                        fit_y = fit_y[:0]
                    if len(fit_x) > 0:
                        _fit_label = 'Lin.Extrap. Fit' if fit_model == 'linear_extrapolated' else 'Linear Fit'
                    ax2.plot(fit_x / time_scale, fit_y, '-', color='red', linewidth=1.5,
                                 label=f'{_fit_label} ({start}, {end_idx})\nR²={fr["R2"]:.3f}')

                    # Dashed red extrapolation line for linear_extrapolated
                    if fr['model'] == 'linear_extrapolated' and len(fit_x) > 0:
                        _t_ro = fr['t_runoff']
                        if np.isfinite(_t_ro) and _t_ro > fit_x[-1]:
                            _extrap_x = np.linspace(fit_x[-1], _t_ro, 50)
                            _extrap_y = _linear_model(_extrap_x, fr['params']['a (slope)'],
                                                      fr['params']['b (intercept)'])
                            ax2.plot(_extrap_x / time_scale, _extrap_y, '--', color='red', linewidth=1.5,
                                     label='_nolegend_')

                    if show_runoff_time and np.isfinite(fr['t_half']) and np.isfinite(fr['t_runoff']):
                        _er = r.get('elong_rate', np.nan)
                        _elong_lbl = f'  ({_er:.2f} aa/s)' if np.isfinite(_er) else ''
                        if fr['model'] != 'linear_extrapolated':
                            ax2.axvline(x=fr['t_half'] / time_scale, color='green', linestyle='--', linewidth=1,
                                        label=fr'$t_{{1/2}}$ ~ {fr["t_half"] / time_scale:.1f} {"h" if time_scale == 60 else "min"}')
                        ax2.axvline(x=fr['t_runoff'] / time_scale, color='orange', linestyle='--', linewidth=1,
                                    label=fr'$\tau_{{runoff}}$ ~ {fr["t_runoff"] / time_scale:.1f} {"h" if time_scale == 60 else "min"}{_elong_lbl}')
                if show_background_line and _bg_raw_value is not None:
                    _bg_y = 0 if remove_background_intensity else _bg_raw_value
                    ax2.axhline(y=_bg_y, color='gray', linestyle=':', linewidth=1,
                                label='_nolegend_')
                ax2.set_xlabel(time_axis_label, fontdict={'size': 16, 'color': 'black'})
                ax2.set_ylabel(y_label, fontdict={'size': 16, 'color': 'black'})
                ax2.tick_params(axis='both', which='major', labelsize=16,
                                labelcolor='black', colors='black')
                for spine in ax2.spines.values():
                    spine.set_color('black')
                    spine.set_linewidth(1.5)
                ax2.set_ylim(ylims)
                if xlims_display is not None:
                    ax2.set_xlim(xlims_display)
                fig2.tight_layout()
                leg2 = ax2.legend(fontsize=10, loc='center left',
                                  bbox_to_anchor=(1.02, 0.5),
                                  framealpha=0.9, edgecolor='black')
                fname = f'{filename_prefix}{plot_name}_ds{ds_idx}_end_{end_idx}'
                fig2.savefig(results_folder / (fname + '.png'), dpi=600,
                             bbox_extra_artists=(leg2,), bbox_inches='tight')
                fig2.savefig(results_folder / (fname + '.svg'), dpi=600,
                             bbox_extra_artists=(leg2,), bbox_inches='tight')
                plt.show()

            print(f'\n── {ds_label} ──')
            _print_fit_table(fit_results, start, fit_model, r2_threshold=r2_threshold,
                             gene_length_effective=gene_length_effective,
                             drug_diffusion_time_min=drug_diffusion_time_min)
            all_fit_results.append(fit_results)

        return all_fit_results

    return fit_results


def _pivot_particle_count_data(particle_counts_df, experiment, value_column,
                               time_values=None):
    """Return time values and a cell-by-time array for one experiment."""
    experiment_df = particle_counts_df[
        particle_counts_df['experiment'] == experiment
    ]
    if experiment_df.empty:
        raise ValueError(f'Experiment {experiment!r} was not found.')
    if value_column not in experiment_df.columns:
        raise ValueError(f'Column {value_column!r} was not found.')
    pivot_df = experiment_df.pivot_table(
        index='cell_key', columns='time_min', values=value_column, aggfunc='first',
    ).sort_index(axis=1)
    if time_values is not None:
        pivot_df = pivot_df.reindex(columns=time_values)
    return pivot_df.columns.to_numpy(dtype=float), pivot_df.to_numpy(dtype=float)


def _prepare_repetition_mean_data(particle_counts_df, repetition_experiments,
                                  value_column, require_all_repetitions):
    """Return time and repetition-mean arrays for one condition."""
    repetition_experiments = tuple(repetition_experiments)
    available_experiments = set(particle_counts_df['experiment'])
    missing_experiments = [
        experiment for experiment in repetition_experiments
        if experiment not in available_experiments
    ]
    if missing_experiments:
        raise ValueError(f'Missing repetitions: {missing_experiments}')

    repetition_df = particle_counts_df[
        particle_counts_df['experiment'].isin(repetition_experiments)
    ]
    repetition_means_df = (
        repetition_df.groupby(['experiment', 'time_min'], sort=False)[value_column]
        .mean()
        .unstack('time_min')
        .reindex(repetition_experiments)
        .sort_index(axis=1)
    )
    if require_all_repetitions:
        repetition_means_df = repetition_means_df.dropna(axis=1, how='any')
    if repetition_means_df.empty or repetition_means_df.shape[1] == 0:
        raise ValueError('No shared repetition time points are available to plot.')
    return (
        repetition_means_df.columns.to_numpy(dtype=float),
        repetition_means_df.to_numpy(dtype=float),
    )


def calculate_actd_cell_half_times(
        particle_counts_df, condition_repetitions,
        value_column='normalized_particle_count', fit_start_idx=0,
        fit_end_idx=None, half_life_target=0.5):
    """Fit each ActD-treated cell and return its target-crossing half-time.

    Each cell is fit independently with ``A * exp(-t / tau) + C`` using the
    same constrained exponential model as the condition-level time-course
    plot. The reported half-time is the duration from the first fitted frame
    to ``half_life_target``; this matches ``half_life_target=0.5`` in
    :func:`plot_multiple_inhibitors`.

    Failed fits and trajectories that cannot reach the requested target are
    retained with a descriptive ``fit_status`` and NaN half-times.
    """
    required_columns = {'condition', 'experiment', 'cell_key', 'time_min', value_column}
    missing_columns = required_columns.difference(particle_counts_df.columns)
    if missing_columns:
        raise ValueError(f'Missing required columns: {sorted(missing_columns)}')
    if not isinstance(condition_repetitions, Mapping) or not condition_repetitions:
        raise ValueError('condition_repetitions must be a non-empty mapping.')

    half_life_target = float(half_life_target)
    if not np.isfinite(half_life_target):
        raise ValueError('half_life_target must be finite.')

    rows = []
    available_experiments = set(particle_counts_df['experiment'])
    for condition, experiments in condition_repetitions.items():
        for experiment in tuple(experiments):
            if experiment not in available_experiments:
                raise ValueError(f'Experiment {experiment!r} was not found.')
            experiment_df = particle_counts_df[
                particle_counts_df['experiment'] == experiment
            ]
            for cell_key, cell_df in experiment_df.groupby('cell_key', sort=False):
                trajectory_df = (
                    cell_df[['time_min', value_column]]
                    .groupby('time_min', as_index=False, sort=True)[value_column]
                    .mean()
                    .sort_values('time_min')
                )
                finite_mask = (
                    np.isfinite(trajectory_df['time_min'].to_numpy(dtype=float))
                    & np.isfinite(trajectory_df[value_column].to_numpy(dtype=float))
                )
                x_data = trajectory_df.loc[finite_mask, 'time_min'].to_numpy(dtype=float)
                y_data = trajectory_df.loc[finite_mask, value_column].to_numpy(dtype=float)
                cell_id = cell_df['cell_id'].iloc[0] if 'cell_id' in cell_df else np.nan

                row = {
                    'condition': condition,
                    'experiment': experiment,
                    'cell_key': cell_key,
                    'cell_id': cell_id,
                    'n_timepoints': len(x_data),
                    'half_life_target': half_life_target,
                    'half_time_min': np.nan,
                    'half_time_h': np.nan,
                    'crossing_time_min': np.nan,
                    'tau_min': np.nan,
                    'plateau': np.nan,
                    'R2': np.nan,
                    'extrapolated': False,
                    'fit_status': 'not_fit',
                }

                start = int(fit_start_idx)
                if len(x_data) < 3:
                    row['fit_status'] = 'insufficient_points'
                    rows.append(row)
                    continue
                if start < 0 or start >= len(x_data):
                    row['fit_status'] = 'invalid_fit_start'
                    rows.append(row)
                    continue

                fit_result = fit_inhibitor_model(
                    x_data,
                    y_data,
                    err_data=None,
                    model='exponential',
                    fit_start_idx=start,
                    fit_end_idx=fit_end_idx,
                )
                if fit_result is None:
                    row['fit_status'] = 'fit_failed'
                    rows.append(row)
                    continue

                params = fit_result['params']
                tau_min = float(params['tau (time constant)'])
                plateau = float(params['C (baseline)'])
                initial_signal = float(y_data[start])
                initial_component = initial_signal - plateau
                target_component = half_life_target - plateau
                row.update({
                    'tau_min': tau_min,
                    'plateau': plateau,
                    'R2': float(fit_result['R2']),
                })

                if initial_component <= 0:
                    row['fit_status'] = 'invalid_initial_signal'
                elif target_component <= 0:
                    row['fit_status'] = 'target_not_reached'
                elif target_component >= initial_component:
                    row['fit_status'] = 'initial_at_or_below_target'
                else:
                    half_time_min = float(
                        tau_min * np.log(initial_component / target_component)
                    )
                    crossing_time_min = float(x_data[start] + half_time_min)
                    observed_duration_min = float(x_data[-1] - x_data[start])
                    row.update({
                        'half_time_min': half_time_min,
                        'half_time_h': half_time_min / 60.0,
                        'crossing_time_min': crossing_time_min,
                        'extrapolated': half_time_min > observed_duration_min,
                        'fit_status': 'valid',
                    })
                rows.append(row)

    return pd.DataFrame.from_records(rows)


def plot_cell_half_times_box_swarm(
        half_times_df, condition_order=None, y_column='half_time_h',
        x_label='', y_label='Half-time (h)', title='', figsize=(6.5, 4.5),
        tick_size=18, swarm_color='black', y_lim=None, show_stats=False,
        calculate_stats=False, max_percentile_significance=99.5,
        x_tick_rotation=0, condition_display_labels=None,
        swarm_max_percentile=None, save_dir=None,
        plot_name='ActD_cell_half_times'):
    """Plot per-cell half-times with the manuscript box-and-swarm style."""
    required_columns = {'condition', y_column}
    missing_columns = required_columns.difference(half_times_df.columns)
    if missing_columns:
        raise ValueError(f'Missing required columns: {sorted(missing_columns)}')

    plot_df = half_times_df.copy()
    if 'fit_status' in plot_df:
        plot_df = plot_df[plot_df['fit_status'] == 'valid']
    plot_df = plot_df[np.isfinite(plot_df[y_column])]
    if plot_df.empty:
        raise ValueError('No valid per-cell half-times are available to plot.')
    if condition_order is None:
        condition_order = list(dict.fromkeys(plot_df['condition']))
    else:
        condition_order = list(condition_order)
    missing_conditions = [
        condition for condition in condition_order
        if not np.any(plot_df['condition'] == condition)
    ]
    if missing_conditions:
        raise ValueError(f'No valid half-times for conditions: {missing_conditions}')

    swarm_df = plot_df
    hidden_swarm_points = 0
    swarm_upper_limit = None
    if swarm_max_percentile is not None:
        swarm_max_percentile = float(swarm_max_percentile)
        if not 0 < swarm_max_percentile <= 100:
            raise ValueError('swarm_max_percentile must be in the interval (0, 100].')
        swarm_upper_limit = float(np.nanpercentile(
            plot_df[y_column].to_numpy(dtype=float), swarm_max_percentile,
        ))
        swarm_df = plot_df[plot_df[y_column] <= swarm_upper_limit]
        hidden_swarm_points = len(plot_df) - len(swarm_df)

    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    sns.boxplot(
        x='condition',
        y=y_column,
        data=plot_df,
        order=condition_order,
        showfliers=False,
        boxprops={'facecolor': 'white', 'edgecolor': 'black'},
        medianprops={'color': 'red'},
        whiskerprops={'color': 'black'},
        capprops={'color': 'black'},
        linewidth=1.5,
        whis=[5, 95],
        width=0.5,
        ax=ax,
    )
    sns.swarmplot(
        x='condition',
        y=y_column,
        data=swarm_df,
        order=condition_order,
        color=swarm_color,
        size=5,
        ax=ax,
    )
    ax.set_facecolor('white')
    ax.set_xlabel(x_label, fontsize=tick_size + 4, fontname='Arial', color='black')
    ax.set_ylabel(y_label, fontsize=tick_size + 4, fontname='Arial', color='black')
    ax.set_title(title, fontsize=tick_size + 4, fontname='Arial', color='black')
    ax.tick_params(axis='x', labelsize=tick_size, colors='black')
    ax.tick_params(axis='y', labelsize=tick_size, colors='black')
    if condition_display_labels is not None:
        ax.set_xticklabels([
            condition_display_labels.get(condition, condition)
            for condition in condition_order
        ])
    for label in ax.get_xticklabels():
        label.set_rotation(x_tick_rotation)
        label.set_fontname('Arial')
        label.set_ha('right' if x_tick_rotation > 0 else 'center')
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')

    global_max = np.nanpercentile(
        plot_df[y_column].to_numpy(dtype=float), max_percentile_significance,
    )
    global_min = np.nanmin(plot_df[y_column].to_numpy(dtype=float))
    global_range = global_max - global_min if global_max != global_min else 1.0
    offset = 0.15 * global_range
    bar_height = 0.02 * global_range
    p_values = {}
    n_bars = 0

    print(f'\n--- Per-cell half-times for {plot_name} ---')
    if hidden_swarm_points:
        print(
            f'  Display: hiding {hidden_swarm_points} swarm point(s) above the '
            f'{swarm_max_percentile:g}th percentile ({swarm_upper_limit:.3f}).'
        )
    for condition in condition_order:
        values = plot_df.loc[plot_df['condition'] == condition, y_column]
        mean_value = float(np.nanmean(values))
        sem_value = float(np.nanstd(values) / np.sqrt(len(values))) if len(values) > 1 else 0.0
        print(f'  {condition}: {mean_value:.3f} ± {sem_value:.3f} (n={len(values)})')
    print('--- End per-cell half-times ---\n')

    if calculate_stats and len(condition_order) > 1:
        for i, condition_1 in enumerate(condition_order[:-1]):
            for j, condition_2 in enumerate(condition_order[i + 1:], start=i + 1):
                group_1 = plot_df.loc[plot_df['condition'] == condition_1, y_column]
                group_2 = plot_df.loc[plot_df['condition'] == condition_2, y_column]
                _, p_value = mannwhitneyu(group_1, group_2)
                comparison_key = f'{condition_1} vs {condition_2}'
                p_values[comparison_key] = float(p_value)
                if p_value < 0.0001:
                    significance = '****'
                elif p_value < 0.001:
                    significance = '***'
                elif p_value < 0.01:
                    significance = '**'
                elif p_value < 0.05:
                    significance = '*'
                else:
                    significance = 'ns'
                if show_stats and significance != 'ns':
                    y_line = global_max + offset * (n_bars + 1)
                    ax.plot(
                        [i, i, j, j],
                        [y_line, y_line + bar_height, y_line + bar_height, y_line],
                        linewidth=1.5,
                        color='black',
                    )
                    ax.text(
                        (i + j) * 0.5,
                        y_line + bar_height,
                        significance,
                        ha='center',
                        va='bottom',
                        color='black',
                        fontsize=tick_size - 4,
                        fontname='Arial',
                    )
                    n_bars += 1

    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_color('black')
    ax.tick_params(
        axis='both', which='major', width=2.0, length=6, colors='black',
    )
    if y_lim is not None and not (show_stats and n_bars > 0):
        ax.set_ylim(y_lim)
    elif show_stats and n_bars > 0:
        highest_bar = global_max + offset * n_bars + bar_height
        padding = 0.1 * (highest_bar - global_min)
        ax.set_ylim(global_min, highest_bar + padding)
    elif swarm_upper_limit is not None:
        lower_limit = 0.0 if global_min >= 0 else global_min
        upper_padding = 0.05 * (swarm_upper_limit - lower_limit)
        ax.set_ylim(lower_limit, swarm_upper_limit + upper_padding)

    fig.tight_layout()
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / f'{plot_name}.png', dpi=600, bbox_inches='tight')
        fig.savefig(save_dir / f'{plot_name}.svg', dpi=600, bbox_inches='tight')
    plt.show()
    return ax, p_values


def plot_actd_conditions(
        particle_counts_df, condition_repetitions, comparison='repetitions',
        control_experiment='control', control_label='Control',
        condition_colors=None, condition_mean_colors=None,
        control_color='#4D4D4D', value_column='normalized_particle_count',
        require_all_repetitions=True, use_sem=True,
        y_label='Normalized Number of Spots', treatment_line_label='5 ug/ml ActD',
        **plot_kwargs):
    """Plot a shared control against one or more ActD condition groups.

    Args:
        particle_counts_df: Output from
            :func:`load_particle_counts_by_experiment`.
        condition_repetitions: Ordered mapping of condition display label to
            its experiment/repetition labels.
        comparison: ``'repetitions'`` plots all repetitions separately;
            ``'average'`` plots one equal-weight mean per condition.
        control_experiment: Experiment label used for the shared control.
        control_label: Control legend label.
        condition_colors: Optional condition-to-repetition-colors mapping.
        condition_mean_colors: Optional condition-to-mean-color mapping.
        control_color: Matplotlib color for the shared control.
        value_column: DataFrame column to plot.
        require_all_repetitions: In average mode, retain only time points
            measured in every repetition within each condition.
        use_sem: If True, show SEM; otherwise show SD.
        y_label: Plot y-axis label.
        treatment_line_label: Label for the optional treatment line.
        **plot_kwargs: Additional options accepted by
            :func:`plot_multiple_inhibitors`.

    Returns:
        Return value from :func:`plot_multiple_inhibitors`.
    """
    comparison = comparison.lower().strip()
    if comparison not in {'repetitions', 'average'}:
        raise ValueError("comparison must be 'repetitions' or 'average'.")
    if not isinstance(condition_repetitions, Mapping) or not condition_repetitions:
        raise ValueError('condition_repetitions must be a non-empty mapping.')

    condition_colors = condition_colors or {}
    condition_mean_colors = condition_mean_colors or {}
    default_palettes = [
        ['#A84300', '#D55E00', '#F3A261'],
        ['#1F4E79', '#3B82B4', '#74A9CF'],
        ['#5B3C88', '#8064A2', '#B39DDB'],
    ]
    control_time, control_data = _pivot_particle_count_data(
        particle_counts_df, control_experiment, value_column,
    )
    full_frames_list = [control_time]
    data_list = [control_data]
    legend_labels = [control_label]
    plot_colors = [control_color]
    # Track total cell counts (for show_n_cells in average mode)
    total_cells_list = [control_data.shape[0]]

    for condition_index, (condition_label, repetitions) in enumerate(
            condition_repetitions.items()):
        repetitions = tuple(repetitions)
        if not repetitions:
            raise ValueError(f'{condition_label!r} has no repetitions.')
        palette = condition_colors.get(
            condition_label,
            default_palettes[condition_index % len(default_palettes)],
        )
        palette = [palette] if isinstance(palette, str) else list(palette)
        if not palette:
            raise ValueError(f'{condition_label!r} has an empty color palette.')

        if comparison == 'repetitions':
            for repetition_index, experiment in enumerate(repetitions, start=1):
                time_values, experiment_data = _pivot_particle_count_data(
                    particle_counts_df, experiment, value_column,
                )
                full_frames_list.append(time_values)
                data_list.append(experiment_data)
                legend_labels.append(f'{condition_label} (Rep. {repetition_index})')
                plot_colors.append(palette[(repetition_index - 1) % len(palette)])
                total_cells_list.append(experiment_data.shape[0])
        else:
            time_values, repetition_mean_data = _prepare_repetition_mean_data(
                particle_counts_df, repetitions, value_column,
                require_all_repetitions=require_all_repetitions,
            )
            full_frames_list.append(time_values)
            data_list.append(repetition_mean_data)
            legend_labels.append(condition_label)
            plot_colors.append(
                condition_mean_colors.get(condition_label, palette[len(palette) // 2])
            )
            # Count total cells across all repetitions for this condition
            _total_cells = sum(
                _pivot_particle_count_data(
                    particle_counts_df, exp, value_column,
                )[1].shape[0]
                for exp in repetitions
            )
            total_cells_list.append(_total_cells)

    plot_kwargs.setdefault('plot_name', f'normalized_spots_conditions_{comparison}')
    plot_kwargs.setdefault('filename_prefix', 'ActD_')
    plot_kwargs.setdefault('show_individual_trajectories', False)
    plot_kwargs.setdefault('show_runoff_time', False)
    plot_kwargs.setdefault('show_treatment_line', False)
    plot_kwargs.setdefault('fit_model', None)
    return plot_multiple_inhibitors(
        full_frames_list=full_frames_list,
        intensities_normalized_list=data_list,
        # inhibitor_frame_index=1: in the ActD spot-counting workflow, each
        # cell is already normalized to its first frame (baseline_frames=1).
        # Frame 0 is the baseline reference; treatment starts at frame 1.
        # This index only affects fitting (fit_start_idx default) and the
        # treatment-line position; both default to off in this function.
        inhibitor_frame_index=1,
        legend_labels=legend_labels,
        colors=plot_colors,
        use_sem=use_sem,
        y_label=y_label,
        treatment_label=treatment_line_label,
        n_cells_override=total_cells_list,
        **plot_kwargs,
    )


def plot_actd_spot_counts(
        particle_counts_df, comparison='repetitions',
        control_experiment='control',
        repetition_experiments=('repetition 1', 'repetition 2', 'repetition 3'),
        control_label='Control', treatment_label='5 ug/ml ActD',
        value_column='normalized_particle_count', require_all_repetitions=True,
        colors=None, use_sem=True, y_label='Normalized Number of Spots',
        **plot_kwargs):
    """Plot a shared control against one ActD condition.

    This backward-compatible wrapper delegates to
    :func:`plot_actd_conditions`.
    """
    comparison = comparison.lower().strip()
    if colors is None:
        colors = (
            ['#4D4D4D', '#D55E00', '#E69F00', '#CC79A7']
            if comparison == 'repetitions'
            else ['#4D4D4D', '#D55E00']
        )
    colors = list(colors)
    if len(colors) < 2:
        raise ValueError('colors must include control and treatment colors.')
    plot_kwargs.setdefault('plot_name', f'normalized_spots_{comparison}')
    return plot_actd_conditions(
        particle_counts_df=particle_counts_df,
        condition_repetitions={treatment_label: tuple(repetition_experiments)},
        comparison=comparison,
        control_experiment=control_experiment,
        control_label=control_label,
        condition_colors={treatment_label: colors[1:]},
        condition_mean_colors={treatment_label: colors[1]},
        control_color=colors[0],
        value_column=value_column,
        require_all_repetitions=require_all_repetitions,
        use_sem=use_sem,
        y_label=y_label,
        treatment_line_label=treatment_label,
        **plot_kwargs,
    )


def process_inhibitor_data(data_dir, inhibitor_frame_index, substring_in_data_dir='', selected_field='spot_int_ch_0', use_sem=True, show_summary=True, max_percentage_threshold_after_treatment=None,
                     frame_interval_sec=60, simulation_dna_sequence=None, inhibitor_delay_time_seconds=60, list_tag_sequences=None, ki_simulation=0.04, ke_simulation=4.5,
                     normalization_method='mean', percentile_range=(5, 95), verbose=False,
                     remove_frame_at_inhibitor_application=False):
    """Process inhibitor runoff experiment data.

    Args:
        data_dir: Path to directory containing results subfolders with
            tracking CSV files.
        inhibitor_frame_index: Frame index at which treatment starts.
        substring_in_data_dir: Filter string to select specific subfolders.
        selected_field: Column name for intensity values
            (e.g., 'spot_int_ch_0').
        use_sem: If True, use SEM for error bars; else use SD.
        show_summary: If True, print summary statistics.
        max_percentage_threshold_after_treatment: Threshold (0-1) to classify
            responding vs non-responding cells.
        frame_interval_sec: Time interval between frames in seconds.
            Default 60 (1 min). Use 20 for 20-second intervals, etc.
        simulation_dna_sequence: DNA sequence for TASEP simulation (optional).
        inhibitor_delay_time_seconds: Delay time for drug to enter cell in
            simulation.
        list_tag_sequences: Tag sequences for probe detection.
        ki_simulation: Initiation rate for simulation.
        ke_simulation: Elongation rate for simulation.
        normalization_method: 'mean' (default) divides by pre-treatment mean
            intensity (values start ~1.0). 'minmax' scales each cell's
            trajectory to [0, 1] range. 'percentile' scales using percentile
            bounds (robust to outliers). None: no normalization.
        percentile_range: Percentile bounds for 'percentile' method.
            Default (5, 95). Use (1, 99) for wider range.
        verbose: If True, print processing details and summary statistics.
            Set to False to suppress all print output.
        remove_frame_at_inhibitor_application: If True, replaces the frame at
            inhibitor application with NaN to remove the focus artifact
            (default False). During live-cell experiments, adding the drug
            mechanically perturbs the sample, causing a transient intensity
            artifact. Setting this to True replaces that frame with NaN.

    Returns:
        Tuple of (responding_indices, time_min_recentered,
        intensities_normalized, array_particles, list_simulation_parameters,
        time_array_sim_min, mean_intensity_ssa_inh, err_intensity_ssa_inh).
    """

    list_dataframes = []
    subfolders = sorted(
        folder for folder in data_dir.iterdir()
        if (folder.is_dir()
            and folder.name.startswith('results_')
            and substring_in_data_dir in folder.name)
    )
    if verbose:
        print('List of processed dataframes:')
    for subfolder in subfolders:
        files = sorted(
            f for f in subfolder.iterdir()
            if (f.is_file()
                and f.name.startswith('tracking_')
                and not f.name.startswith('._')
                and f.suffix.lower() == '.csv')
        )
        if len(files) != 1:
            raise ValueError(
                f'{subfolder.name} must contain exactly one real tracking_*.csv '
                f'file after excluding AppleDouble sidecars; found {len(files)}.'
            )
        dfs = pd.read_csv(files[0])
        list_dataframes.append(dfs)
        if verbose:
            print('subfolder:', subfolder)

    # detect the maximum frame number across all dataframes
    max_frame = 0
    for df in list_dataframes:
        max_frame = max(max_frame, df['frame'].max())

    # terminate the program if the maximum frame is less than 1
    if max_frame < 10:
        if verbose:
            print('No dataframes found with frame number greater than 10 frames.')
        return None, None, None, None, None, None, None, None

    # frame_indices includes max_frame so that reindex() covers every
    # observed frame (0 through max_frame inclusive).  The returned
    # time_min_recentered and intensities_normalized are both derived from
    # frame_indices, keeping them aligned.
    frame_indices = np.arange(0, max_frame + 1)
    array_particles = np.zeros((len(list_dataframes), len(frame_indices)))
    average_intensity = np.zeros((len(list_dataframes), len(frame_indices)))
    intensities_normalized = np.zeros((len(list_dataframes), len(frame_indices)))
    average_molecules_before_treatment = np.zeros(len(list_dataframes))
    for i, df in enumerate(list_dataframes):
        # Group by 'frame' and count unique particles
        particle_counts_per_frame = df.groupby('frame')['particle'].nunique()
        sum_intensities_per_frame = df.groupby('frame')[selected_field].sum()
        # Re-index to include frames with no particles (fill missing with 0)
        particle_counts_per_frame = particle_counts_per_frame.reindex(frame_indices, fill_value=0).values
        sum_intensities_per_frame = sum_intensities_per_frame.reindex(frame_indices, fill_value=0).values
        # Ensure the lengths match frame_indices
        if len(particle_counts_per_frame) > len(frame_indices):
            particle_counts_per_frame = particle_counts_per_frame[:len(frame_indices)]
            sum_intensities_per_frame = sum_intensities_per_frame[:len(frame_indices)]
        # Compute raw intensities (always use 'mean' per-cell first)
        intensities_normalized_before_treatment_intensity, average_intensity_with_respect_number_particles, average_particles_before_treatment = calculate_intensity(
            particle_counts_per_frame, sum_intensities_per_frame, inhibitor_frame_index,
            normalization_method=normalization_method
        )
        # Store in array
        normalized_particles, _ = calculate_number_of_particles_per_frame(particle_counts_per_frame, inhibitor_frame_index)
        array_particles[i] = normalized_particles
        average_intensity[i] = average_intensity_with_respect_number_particles
        intensities_normalized[i] = intensities_normalized_before_treatment_intensity
        average_molecules_before_treatment[i] = average_particles_before_treatment

    # Apply global normalization across all cells (after the loop)
    if normalization_method == 'minmax':
        global_min = intensities_normalized.min()
        global_max = intensities_normalized.max()
        if global_max - global_min > 0:
            intensities_normalized = (intensities_normalized - global_min) / (global_max - global_min)
        if verbose:
            print(f'Global minmax normalization applied (min={global_min:.4f}, max={global_max:.4f})')
    elif normalization_method == 'percentile':
        global_low = np.percentile(intensities_normalized, percentile_range[0])
        global_high = np.percentile(intensities_normalized, percentile_range[1])
        if global_high - global_low > 0:
            intensities_normalized = (intensities_normalized - global_low) / (global_high - global_low)
        if verbose:
            print(f'Global percentile normalization applied (P{percentile_range[0]}={global_low:.4f}, P{percentile_range[1]}={global_high:.4f})')

    # ── Remove artifact frame at inhibitor application ────────────────
    # During live-cell experiments, adding the inhibitor (e.g., pipetting
    # harringtonine into the dish) mechanically perturbs the sample,
    # causing cells to transiently go out of focus. This creates an
    # artificial intensity dip at the treatment frame that does not
    # reflect actual translational run-off. Replacing this frame with
    # NaN removes the artifact while preserving the time axis, so the
    # mean, SEM, and model fits are not biased by the focus disturbance.
    if remove_frame_at_inhibitor_application:
        intensities_normalized[:, inhibitor_frame_index] = np.nan
        array_particles[:, inhibitor_frame_index] = np.nan
        if verbose:
            print(f'Removed frame at inhibitor application '
                  f'(index {inhibitor_frame_index}) → replaced with NaN')

    treatment_start_index = inhibitor_frame_index
    non_responding_indices = []  # average post-treatment is not decreasing below the threshold
    responding_indices = []
    # Classify each trajectory
    if max_percentage_threshold_after_treatment is not None:
        for i, trajectory in enumerate(intensities_normalized):
            baseline = np.nanmean(trajectory[:treatment_start_index])
            threshold = max_percentage_threshold_after_treatment * baseline
            avg_post_treatment = np.nanmean(trajectory[treatment_start_index:])
            if avg_post_treatment >= threshold:
                non_responding_indices.append(i)
            else:
                responding_indices.append(i)
    else:
        # print a warning if the threshold is not provided
        if verbose:
            print('Warning: No threshold provided. All trajectories are considered responding.')
        # if no threshold is provided, all trajectories are considered responding
        responding_indices = list(range(len(intensities_normalized)))
        non_responding_indices = []

    if show_summary and verbose:
        print(f'\n{"─" * 40}')
        print('Summary')
        print(f'{"─" * 40}')
        print('\nCells')
        print(f'{"─" * 40}')
        # report total cells
        total_cells = len(intensities_normalized)
        print('total_cells: ', total_cells)
        number_of_responding_cells = len(responding_indices)
        print('number_of_responding_cells: ', number_of_responding_cells)
        number_of_non_responding_cells = len(non_responding_indices)
        print('number_of_non_responding_cells: ', number_of_non_responding_cells)
        # Denominator equals total_cells (every trajectory is classified as
        # responding or non-responding), which is guaranteed > 0 here because
        # the function returns early when max_frame < 10.
        percentage_non_responding_cells = number_of_non_responding_cells / total_cells * 100
        print('percentage_non_responding_cells: ', np.round(percentage_non_responding_cells,1))
        print(f'\nRNA\n{"─" * 40}')
        total_molecules_before_treatment = average_molecules_before_treatment.sum()
        print('average_total_molecules_before_treatment: ', total_molecules_before_treatment)
        average_molecules_before_treatment_responding = average_molecules_before_treatment[responding_indices].sum()
        print('average_molecules_before_treatment_responding: ', np.round(average_molecules_before_treatment_responding,1))
        average_molecules_before_treatment_non_responding = average_molecules_before_treatment[non_responding_indices].sum()
        print('average_molecules_before_treatment_non_responding: ', average_molecules_before_treatment_non_responding)

    # Convert frame indices to time in minutes, centered on treatment
    time_min_recentered = (frame_indices - inhibitor_frame_index) * frame_interval_sec / 60.0


    if simulation_dna_sequence is not None:
        if list_tag_sequences is None:
            list_tag_sequences = [HA_TAG]
        list_simulation_parameters = [ki_simulation, ke_simulation]
        number_repetitions = 100
        burnin_time = 2000
        t_max = max_frame * frame_interval_sec  # Total experiment time in seconds
        step_size_in_sec = 1
        time_array = np.arange(0, t_max, step_size_in_sec)
        _, rna, gene_length, first_probe_position_vector, second_probe_position_vector = read_gene_sequence_return_probes(simulation_dna_sequence, min_protein_length=50, list_tag_sequences=list_tag_sequences)
        ke = calculate_codon_elongation_rates(rna, global_elongation_rate=list_simulation_parameters[1])
        time_perturbation_application = inhibitor_frame_index * frame_interval_sec + inhibitor_delay_time_seconds
        evaluating_inhibitor = 1
        # inhibitor_effectiveness=95: percentage reduction in initiation rate
        # upon drug application.  The 95% value is specific to the
        # translational-inhibitor run-off assay modeled here (see Tanenbaum
        # et al. 2015).  The optimization/Harringtonine validation uses 90%
        # in a different experimental context.
        ssa_array = simulate_TASEP_SSA(list_simulation_parameters[0],
                                    ke,
                                    gene_length,
                                    t_max,
                                    time_interval_in_seconds=step_size_in_sec,
                                    number_repetitions=number_repetitions,
                                    first_probe_position_vector=first_probe_position_vector,
                                    timePerturbationApplication=time_perturbation_application,
                                    inhibitor_effectiveness=95,
                                    evaluatingInhibitor=evaluating_inhibitor,
                                    burnin_time=burnin_time,
                                    constant_elongation_rate=list_simulation_parameters[1],
                                    fast_output=True)[2]
        # downsample time array and simulation output to match experimental frame interval
        downsample_factor = int(frame_interval_sec)  # e.g., 60 for 1 min, 20 for 20 sec
        time_array_downsampled = time_array[::downsample_factor]
        ssa_array_downsampled = ssa_array[:, ::downsample_factor]
        # plotting
        time_array_sim_min = time_array_downsampled/60
        normalized_data = np.zeros_like(ssa_array_downsampled)
        for i in range(ssa_array_downsampled.shape[0]):
            mean_before_treatment = np.mean(ssa_array_downsampled[i,:inhibitor_frame_index])
            if mean_before_treatment == 0:
                normalized_data[i] = np.zeros_like(ssa_array_downsampled[i])
            else:
                normalized_data[i] = ssa_array_downsampled[i]/mean_before_treatment
        mean_intensity_ssa_inh = np.mean(normalized_data, axis=0)
        if use_sem:
            err_intensity_ssa_inh = np.std(normalized_data, axis=0) / np.sqrt(normalized_data.shape[0])
        else:
            err_intensity_ssa_inh = np.std(normalized_data, axis=0)
    else:
        mean_intensity_ssa_inh = None
        err_intensity_ssa_inh = None
        time_array_sim_min = None
        list_simulation_parameters = [None, None]

    return responding_indices, time_min_recentered, intensities_normalized, array_particles, list_simulation_parameters, time_array_sim_min, mean_intensity_ssa_inh, err_intensity_ssa_inh


def simulate_inhibitor(gene_sequence, ki=0.04, ke_global=5, use_sem=False, max_frame=20,
                       list_tag_sequences=None, inhibitor_frame=5, inhibitor_delay_seconds=60,
                       verbose=False):
    """Run a TASEP simulation for a gene sequence with inhibitor treatment.

    Simulation-only function (no experimental data). Uses codon-usage-aware
    elongation rates and supports multi-tag probe detection.

    Args:
        gene_sequence: DNA sequence to simulate.
        ki: Initiation rate.
        ke_global: Global elongation rate.
        use_sem: If True, error bars use SEM; else SD.
        max_frame: Maximum number of frames (minutes) to simulate.
        list_tag_sequences: Tag sequences for probe detection.
            Defaults to [HA_TAG, GFP_TAG].
        inhibitor_frame: Frame (minute) at which inhibitor treatment starts.
        inhibitor_delay_seconds: Delay time for drug to enter cell in
            simulation (seconds).
        verbose: If True, print probe positions and other diagnostics.


    Returns:
        Tuple of (time_array_min, mean_intensity_ssa_inh,
        err_intensity_ssa_inh).
    """
    if list_tag_sequences is None:
        list_tag_sequences = [HA_TAG, GFP_TAG]

    # reading the gene sequence
    protein, rna, _, indexes_tags, _, _ ,_ = read_sequence(seq=gene_sequence, min_protein_length=50, TAG=list_tag_sequences)
    gene_length = len(protein)
    tag_positions_first_probe_vector = indexes_tags[0]
    first_probe_position_vector = create_probe_vector(tag_positions_first_probe_vector, gene_length)

    # second probe vector.
    tag_positions_second_probe_vector = indexes_tags[1]
    second_probe_position_vector = create_probe_vector(tag_positions_second_probe_vector, gene_length)

    if verbose:
        print(f"First probe positions: {tag_positions_first_probe_vector}")
        print(f"Second probe positions: {tag_positions_second_probe_vector}")

    ke_codon_dependent = calculate_codon_elongation_rates(rna, global_elongation_rate=ke_global)

    full_frames = np.arange(0, max_frame)
    full_frames = full_frames - inhibitor_frame

    # This function assumes 1-minute frame intervals by design (see docstring:
    # "max_frame: Maximum number of frames (minutes)").  For non-60-second
    # intervals, use process_inhibitor_data which accepts frame_interval_sec.
    number_repetitions = 100
    burnin_time = 2000
    t_max = max_frame * 60
    step_size_in_sec = 1
    time_array = np.arange(0, t_max, step_size_in_sec)
    time_perturbation_application = inhibitor_frame * 60 + inhibitor_delay_seconds
    evaluating_inhibitor = 1
    ssa_array = simulate_TASEP_SSA(ki,
                                ke_codon_dependent,
                                gene_length,
                                t_max,
                                time_interval_in_seconds=step_size_in_sec,
                                number_repetitions=number_repetitions,
                                first_probe_position_vector=first_probe_position_vector,
                                timePerturbationApplication=time_perturbation_application,
                                inhibitor_effectiveness=95,
                                evaluatingInhibitor=evaluating_inhibitor,
                                burnin_time=burnin_time,
                                constant_elongation_rate=None,  # codon-usage aware
                                fast_output=False)[2]
    # downsample time array and simulation output to 1 minute resolution.
    time_array_downsampled = time_array[::60]
    ssa_array_downsampled = ssa_array[:, ::60]
    # plotting
    time_array_min = time_array_downsampled / 60
    normalized_data = np.zeros_like(ssa_array_downsampled)
    for i in range(ssa_array_downsampled.shape[0]):
        mean_before_treatment = np.mean(ssa_array_downsampled[i, :inhibitor_frame])
        if mean_before_treatment == 0:
            normalized_data[i] = np.zeros_like(ssa_array_downsampled[i])
        else:
            normalized_data[i] = ssa_array_downsampled[i] / mean_before_treatment
    mean_intensity_ssa_inh = np.mean(normalized_data, axis=0)

    if use_sem:
        err_intensity_ssa_inh = np.std(normalized_data, axis=0) / np.sqrt(normalized_data.shape[0])
    else:
        err_intensity_ssa_inh = np.std(normalized_data, axis=0)

    return time_array_min, mean_intensity_ssa_inh, err_intensity_ssa_inh


def plot_inhibitor_simulation(legend_label, time_array_min, mean_intensity_ssa_inh, err_intensity_ssa_inh,
                              figsize=(6, 3), ylims=(0, 1.4), inhibitor_frame=5):
    """Plot a single inhibitor simulation result.

    Args:
        legend_label: Label for the simulation trace.
        time_array_min: Time array in minutes.
        mean_intensity_ssa_inh: Mean normalized intensity from simulation.
        err_intensity_ssa_inh: Error (SD or SEM) of normalized intensity.
        figsize: Figure size.
        ylims: Y-axis limits.
        inhibitor_frame: Frame (minute) at which inhibitor treatment starts
            (for time offset).
    """
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    ax.set_facecolor('white')

    # Plotting the model
    ax.plot(time_array_min - inhibitor_frame, mean_intensity_ssa_inh, '-', color='red', linewidth=4, label=legend_label)
    ax.fill_between(time_array_min - inhibitor_frame, mean_intensity_ssa_inh - err_intensity_ssa_inh,
                    mean_intensity_ssa_inh + err_intensity_ssa_inh, color='red', alpha=0.1)

    # Plot the inhibitor frame as a vertical line at x=0
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1)

    ax.set_xlabel("Time (min)", fontdict={'size': 16})
    ax.set_ylabel("Norm. Intensity", fontdict={'size': 16})
    ax.tick_params(axis='both', which='major', labelsize=12)

    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)

    ax.set_ylim(ylims)
    fig.tight_layout()
    ax.legend(fontsize=10)
    plt.show()


def plot_multiple_inhibitor_simulations(
    legend_labels,
    time_arrays_min,
    mean_intensities,
    err_intensities,
    figsize=(8, 5),
    ylims=(0, 1.5),
    inhibitor_frame=5
):
    """Overlay multiple inhibitor simulation results on a single plot.

    Args:
        legend_labels: Labels for each simulation trace (str or list of str).
        time_arrays_min: Time arrays in minutes (array or list of arrays).
        mean_intensities: Mean normalized intensities from each simulation.
        err_intensities: Errors (SD or SEM) for each simulation.
        figsize: Figure size.
        ylims: Y-axis limits.
        inhibitor_frame: Frame (minute) at which inhibitor treatment starts
            (for time offset).
    """
    # --- normalize legend_labels to list ---
    if not isinstance(legend_labels, (list, tuple)):
        legend_labels = [legend_labels]
    n = len(legend_labels)

    # --- broadcast time_arrays_min if needed ---
    if not isinstance(time_arrays_min, (list, tuple)):
        time_arrays_min = [time_arrays_min] * n

    # --- similarly ensure mean_intensities and err_intensities are lists ---
    if not isinstance(mean_intensities, (list, tuple)):
        mean_intensities = [mean_intensities] * n
    if not isinstance(err_intensities, (list, tuple)):
        err_intensities = [err_intensities] * n

    # sanity check
    if not (len(time_arrays_min) == n == len(mean_intensities) == len(err_intensities)):
        raise ValueError('All four inputs must be lists of the same length.')

    # 1) create figure
    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    ax.set_facecolor('white')

    # 2) choose a colormap and generate n distinct colors
    cmap = plt.get_cmap('tab10')
    colors = cmap(np.linspace(0, 1, n))

    # 3) plot each simulation with its assigned color
    for idx, (label, t_min, mean_inh, err_inh) in enumerate(zip(
        legend_labels, time_arrays_min, mean_intensities, err_intensities
    )):
        t_plot = np.array(t_min) - inhibitor_frame
        color = colors[idx]
        ax.plot(t_plot, mean_inh,
                '-', color=color, linewidth=3, label=label)
        ax.fill_between(
            t_plot,
            mean_inh - err_inh,
            mean_inh + err_inh,
            color=color,
            alpha=0.05
        )

    # vertical line at t=0
    ax.axvline(x=0, linestyle='--', color='k', linewidth=1)

    # axes labels and formatting
    ax.set_xlabel("Time (min)", size=16)
    ax.set_ylabel("Norm. Intensity", size=16)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
    ax.set_ylim(ylims)

    # legend and layout
    ax.legend(fontsize=10, frameon=False)
    fig.tight_layout()
    plt.show()
