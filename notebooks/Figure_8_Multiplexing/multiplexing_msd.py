"""Paired-channel MSD analysis used by Multiplexing_MSD_Analysis.ipynb.

This path-agnostic module loads validated MicroLive tracking exports, computes
cell-level normal and anomalous MSD fits, summarizes matched-channel spot
properties, performs paired channel statistics, and creates the notebook's
publication figures. Source files are pooled; cells are the biological samples
and tracking channels are paired conditions within each cell.
"""

from __future__ import annotations

import math
import re
import warnings
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import trackpy as tp
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from microlive import microscopy as mi
from scipy.stats import linregress, ttest_rel, wilcoxon

_BASE_TABLE_FILENAMES = {
    "dataset_inventory": "dataset_inventory.csv",
    "trajectory_filter_summary": "trajectory_filter_summary.csv",
    "trajectory_summary": "trajectory_summary.csv",
    "particle_msd": "particle_msd_long.csv",
    "cell_msd": "cell_msd_curves.csv",
    "fit_summary": "msd_fit_summary.csv",
    "cell_summary": "cell_channel_summary.csv",
}

_CHANNEL_COLORS = (
    "#00CC33",  # green (channel 0)
    "#FF00FF",  # magenta (channel 1)
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
)
_CHANNEL_MARKERS = ("o", "s", "^", "D")

_FINAL_TABLE_FILENAMES = {
    "trajectory_spot_summary": "trajectory_spot_summary.csv",
    "cell_spot_summary": "cell_spot_summary.csv",
    "cell_analysis_summary": "cell_analysis_summary.csv",
    "channel_condition_summary": "channel_condition_summary.csv",
    "msd_d_consistency": "msd_d_consistency.csv",
    "paired_channel_statistics": "paired_channel_statistics.csv",
}

METRIC_SPECS = {
    "D_um2_s": {
        "label": "Diffusion coefficient, D (µm²/s)",
        "title": "Normal diffusion",
        "unit": "µm²/s",
        "trajectory_metric": None,
        "default_log": False,
        "zero_baseline": True,
    },
    "alpha": {
        "label": "Anomalous exponent, α",
        "title": "Anomalous diffusion exponent",
        "unit": "dimensionless",
        "trajectory_metric": None,
        "default_log": False,
        "zero_baseline": True,
    },
    "K_alpha_um2_s_alpha": {
        "label": "Generalized diffusion coefficient, Kα (µm²/s^α)",
        "title": "Anomalous generalized diffusion",
        "unit": "µm²/s^α",
        "trajectory_metric": None,
        "default_log": True,
        "zero_baseline": False,
    },
    "cell_spot_intensity_au": {
        "label": "Mean spot intensity (a.u.)",
        "title": "Spot intensity",
        "unit": "a.u.",
        "trajectory_metric": "trajectory_spot_intensity_au",
        "default_log": False,
        "zero_baseline": False,
    },
    "cell_spot_snr": {
        "label": "Mean spot SNR",
        "title": "Spot signal-to-noise ratio",
        "unit": "dimensionless",
        "trajectory_metric": "trajectory_spot_snr",
        "default_log": False,
        "zero_baseline": True,
    },
    "cell_spot_fwhm_um": {
        "label": "Mean spot FWHM (µm)",
        "title": "Spot size",
        "unit": "µm",
        "trajectory_metric": "trajectory_spot_fwhm_um",
        "default_log": False,
        "zero_baseline": False,
    },
    "cell_spot_fwhm_nm": {
        "label": "Mean spot FWHM (nm)",
        "title": "Spot size",
        "unit": "nm",
        "trajectory_metric": "trajectory_spot_fwhm_nm",
        "default_log": False,
        "zero_baseline": False,
    },
    "n_trajectories": {
        "label": "No. trajectories per cell",
        "title": "Tracking yield",
        "unit": "trajectories",
        "trajectory_metric": None,
        "default_log": False,
        "zero_baseline": True,
    },
    "median_trajectory_frames": {
        "label": "Median trajectory length (frames)",
        "title": "Trajectory duration",
        "unit": "frames",
        "trajectory_metric": None,
        "default_log": False,
        "zero_baseline": True,
    },
}


# Data loading and MSD fitting


def _metadata_value(text: str, field_name: str, metadata_path: Path) -> str:
    """Extract a required metadata field from a text file."""
    pattern = rf"^\s*{re.escape(field_name)}\.+\s*(.*?)\s*$"
    match = re.search(pattern, text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"{field_name!r} was not found in {metadata_path}.")
    return match.group(1).strip()


def _optional_metadata_value(
    text: str,
    field_name: str,
    metadata_path: Path,
    default: str | None = None,
) -> str | None:
    """Extract an optional metadata field, returning a default when absent."""
    try:
        return _metadata_value(text, field_name, metadata_path)
    except ValueError:
        return default


def _parse_channel_values(value: str | None) -> tuple[int, ...]:
    """Parse unique integer channel identifiers from metadata text."""
    if not value:
        return ()
    return tuple(dict.fromkeys(int(item) for item in re.findall(r"-?\d+", value)))


def _read_microlive_metadata(metadata_path: str | Path) -> dict[str, Any]:
    """Read the calibration and acquisition fields needed for MSD analysis."""
    path = Path(metadata_path).expanduser().resolve()
    text = path.read_text(encoding="utf-8", errors="replace")

    frame_interval_s = float(_metadata_value(text, "Time Interval (s)", path))
    active_frame_count = int(_metadata_value(text, "Active Frame Count", path))
    xy_nm_per_px = float(_metadata_value(text, "Voxel Size YX (nm)", path))
    z_text = _optional_metadata_value(text, "Voxel Size Z (nm)", path)
    z_nm_per_px = float(z_text) if z_text not in (None, "", "None") else np.nan
    projection_mode = _optional_metadata_value(text, "Projection Mode", path, "Unknown")
    tracked_channels = _parse_channel_values(
        _optional_metadata_value(text, "Tracked Channels", path)
    )
    snr_method_by_channel = {
        int(channel): method.strip()
        for channel, method in re.findall(
            r"^\s*Channel\s+(\d+)\s+SNR Method\.+\s*(.*?)\s*$",
            text,
            flags=re.MULTILINE,
        )
    }

    if not np.isfinite(frame_interval_s) or frame_interval_s <= 0:
        raise ValueError(f"Invalid frame interval in {path}.")
    if active_frame_count <= 1:
        raise ValueError(f"Active frame count must be at least two in {path}.")
    if not np.isfinite(xy_nm_per_px) or xy_nm_per_px <= 0:
        raise ValueError(f"Invalid XY voxel size in {path}.")

    projection_lower = str(projection_mode).lower()
    is_3d = "3d" in projection_lower and "2d" not in projection_lower
    if is_3d and (not np.isfinite(z_nm_per_px) or z_nm_per_px <= 0):
        raise ValueError(
            f"True 3D tracking requires a positive Z voxel size in {path}."
        )

    return {
        "metadata_path": str(path),
        "source_lif_path": _optional_metadata_value(text, "Data Folder Path", path, ""),
        "selected_image": _optional_metadata_value(text, "Image Name", path, ""),
        "frame_interval_s": frame_interval_s,
        "active_frame_count": active_frame_count,
        "xy_um_per_px": xy_nm_per_px / 1000.0,
        "z_um_per_px": z_nm_per_px / 1000.0 if np.isfinite(z_nm_per_px) else np.nan,
        "projection_mode": projection_mode,
        "tracked_channels_metadata": tracked_channels,
        "snr_method_by_channel": snr_method_by_channel,
        "is_3d": is_3d,
    }


def _infer_condition(source_id: str) -> str:
    """Infer an experimental condition from a result-folder identifier."""
    source_upper = source_id.upper()
    for condition in ("DISH2", "KDM5B"):
        if condition in source_upper:
            return condition
    cleaned = re.sub(r"^results_", "", source_id, flags=re.IGNORECASE)
    cleaned = re.sub(r"^\d{8}_", "", cleaned)
    cleaned = re.split(r"_Position|_Series|_multiplexing", cleaned, maxsplit=1)[0]
    return cleaned or "Unassigned"


def _validate_time_column(
    dataframe: pd.DataFrame,
    frame_interval_s: float,
    csv_path: Path,
) -> None:
    """Validate tracking time values against the metadata frame interval."""
    if "time" not in dataframe.columns or len(dataframe) < 2:
        return
    frame_time = dataframe[["frame", "time"]].apply(pd.to_numeric, errors="coerce")
    frame_time = frame_time.replace([np.inf, -np.inf], np.nan).dropna()
    if frame_time.empty:
        raise ValueError(f"{csv_path.name} has no finite tracking times.")
    per_frame = (
        frame_time.groupby("frame", as_index=False)["time"]
        .median()
        .sort_values("frame")
    )
    frame_delta = np.diff(per_frame["frame"].to_numpy(dtype=float))
    time_delta = np.diff(per_frame["time"].to_numpy(dtype=float))
    valid = frame_delta > 0
    if not np.any(valid):
        return
    observed = float(np.median(time_delta[valid] / frame_delta[valid]))
    tolerance = max(1e-9, abs(frame_interval_s) * 1e-6)
    if not np.isclose(observed, frame_interval_s, rtol=1e-6, atol=tolerance):
        raise ValueError(
            f"{csv_path.name} reports {observed:g} s/frame in its time column, "
            f"but its metadata reports {frame_interval_s:g} s/frame."
        )


def _normalize_tracking_dataframe(
    dataframe: pd.DataFrame,
    *,
    source_id: str,
    condition: str,
    csv_path: Path,
    metadata: Mapping[str, Any],
) -> pd.DataFrame:
    """Validate and standardize one tracking dataframe."""
    required = {"cell_id", "spot_type", "frame", "x", "y"}
    missing = sorted(required.difference(dataframe.columns))
    if (
        "particle" not in dataframe.columns
        and "unique_particle" not in dataframe.columns
    ):
        missing.append("particle or unique_particle")
    if metadata["is_3d"] and "z" not in dataframe.columns:
        missing.append("z")
    if missing:
        raise ValueError(f"{csv_path.name} is missing columns: {', '.join(missing)}.")
    if dataframe.empty:
        raise ValueError(f"{csv_path.name} is empty.")
    if dataframe["cell_id"].isna().any() or dataframe["spot_type"].isna().any():
        raise ValueError(
            f"{csv_path.name} contains missing cell or channel identifiers."
        )

    normalized = dataframe.copy()
    numeric_columns = ["frame", "x", "y", "spot_type"]
    if metadata["is_3d"]:
        numeric_columns.append("z")
    numeric = normalized[numeric_columns].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"{csv_path.name} contains non-finite coordinates or frames.")
    if not np.allclose(
        numeric["frame"].to_numpy(dtype=float),
        np.round(numeric["frame"].to_numpy(dtype=float)),
        rtol=0,
        atol=1e-9,
    ):
        raise ValueError(f"{csv_path.name} contains non-integer frame indices.")
    if not np.allclose(
        numeric["spot_type"].to_numpy(dtype=float),
        np.round(numeric["spot_type"].to_numpy(dtype=float)),
        rtol=0,
        atol=1e-9,
    ):
        raise ValueError(f"{csv_path.name} contains non-integer channel identifiers.")

    normalized["frame"] = numeric["frame"].astype(int)
    normalized["spot_type"] = numeric["spot_type"].astype(int)
    normalized["x"] = numeric["x"]
    normalized["y"] = numeric["y"]
    if metadata["is_3d"]:
        normalized["z"] = numeric["z"]
    if (normalized["frame"] < 0).any():
        raise ValueError(f"{csv_path.name} contains negative frame indices.")
    if (normalized["frame"] >= int(metadata["active_frame_count"])).any():
        raise ValueError(
            f"{csv_path.name} contains frames outside its metadata frame range."
        )

    _validate_time_column(normalized, float(metadata["frame_interval_s"]), csv_path)

    particle_column = (
        "particle" if "particle" in normalized.columns else "unique_particle"
    )
    if normalized[particle_column].isna().any():
        raise ValueError(f"{csv_path.name} contains missing particle identifiers.")

    normalized["condition"] = condition
    normalized["source_id"] = source_id
    normalized["local_cell_id"] = normalized["cell_id"].astype(str)
    normalized["channel"] = normalized["spot_type"].astype(int)
    normalized["original_particle_id"] = normalized[particle_column].astype(str)
    normalized["global_cell_id"] = (
        normalized["source_id"] + "::cell_" + normalized["local_cell_id"]
    )
    normalized["global_particle_id"] = (
        normalized["global_cell_id"]
        + "::ch_"
        + normalized["channel"].astype(str)
        + "::particle_"
        + normalized["original_particle_id"]
    )

    duplicate_mask = normalized.duplicated(["global_particle_id", "frame"], keep=False)
    if duplicate_mask.any():
        example = normalized.loc[duplicate_mask, ["global_particle_id", "frame"]].iloc[
            0
        ]
        raise ValueError(
            f"{csv_path.name} has duplicate detections for "
            f"{example['global_particle_id']} at frame {example['frame']}."
        )
    return normalized


def _load_results(
    results_root: str | Path,
    *,
    strict: bool = True,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and validate all ``results_*`` folders without modifying source data."""
    root = Path(results_root).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Results directory not found: {root}")
    result_folders = sorted(
        folder
        for folder in root.iterdir()
        if folder.is_dir() and folder.name.startswith("results_")
    )
    if not result_folders:
        raise ValueError(f"No results_* folders found in {root}.")

    inventory_rows: list[dict[str, Any]] = []
    tracking_parts: list[pd.DataFrame] = []
    for folder in result_folders:
        tracking_files = sorted(
            path
            for path in folder.glob("tracking_*.csv")
            if not path.name.startswith("._")
        )
        metadata_files = sorted(
            path
            for path in folder.glob("Metadata_*.txt")
            if not path.name.startswith("._")
        )
        row: dict[str, Any] = {
            "condition": _infer_condition(folder.name),
            "source_id": folder.name,
            "result_folder": str(folder),
            "tracking_csv": str(tracking_files[0]) if len(tracking_files) == 1 else "",
            "metadata_path": str(metadata_files[0]) if len(metadata_files) == 1 else "",
            "validation_status": "invalid",
            "validation_reason": "",
        }
        try:
            if len(tracking_files) != 1 or len(metadata_files) != 1:
                raise ValueError(
                    f"Expected one tracking_*.csv and one Metadata_*.txt; found "
                    f"{len(tracking_files)} and {len(metadata_files)}."
                )
            metadata = _read_microlive_metadata(metadata_files[0])
            raw = pd.read_csv(tracking_files[0])
            normalized = _normalize_tracking_dataframe(
                raw,
                source_id=folder.name,
                condition=row["condition"],
                csv_path=tracking_files[0],
                metadata=metadata,
            )
            for key, value in metadata.items():
                if key not in {"metadata_path", "snr_method_by_channel"}:
                    normalized[key] = [value] * len(normalized)
            normalized["snr_method"] = normalized["channel"].map(
                metadata["snr_method_by_channel"]
            )
            channels = tuple(sorted(normalized["channel"].unique().tolist()))
            row.update(
                {
                    **metadata,
                    "n_rows": len(normalized),
                    "n_local_cells": int(normalized["global_cell_id"].nunique()),
                    "channels_csv": ",".join(map(str, channels)),
                    "first_frame": int(normalized["frame"].min()),
                    "last_frame": int(normalized["frame"].max()),
                    "validation_status": "valid",
                    "validation_reason": "",
                }
            )
            tracking_parts.append(normalized)
        except Exception as exc:
            row["validation_reason"] = str(exc)
        inventory_rows.append(row)

    inventory = pd.DataFrame(inventory_rows)
    invalid_rows = [
        row for row in inventory_rows if row["validation_status"] != "valid"
    ]
    if invalid_rows and strict:
        reasons = "\n".join(
            f"- {row['source_id']}: {row['validation_reason']}"
            for row in invalid_rows
        )
        raise ValueError(f"Invalid tracking result folders:\n{reasons}")
    if invalid_rows:
        reasons = "; ".join(
            f"{row['source_id']}: {row['validation_reason']}"
            for row in invalid_rows
        )
        warnings.warn(
            f"Continuing without {len(invalid_rows)} invalid result folder(s): "
            f"{reasons}",
            RuntimeWarning,
            stacklevel=2,
        )
    if not tracking_parts:
        reasons = "\n".join(
            f"- {row['source_id']}: {row['validation_reason']}"
            for row in inventory_rows
        )
        raise ValueError(f"No valid tracking datasets were loaded:\n{reasons}")
    tracking = pd.concat(tracking_parts, ignore_index=True, sort=False)

    if verbose:
        valid_count = int((inventory["validation_status"] == "valid").sum())
        invalid_count = int((inventory["validation_status"] != "valid").sum())
        channels = sorted(tracking["channel"].unique().tolist())
        print(
            f"Loaded {valid_count} valid result folders "
            f"({invalid_count} invalid); channels={channels}."
        )
    return inventory, tracking


def _trajectory_summary_for_unit(unit: pd.DataFrame) -> pd.DataFrame:
    """Summarize trajectory lengths and observation windows for one unit."""
    metadata = unit.iloc[0]
    dt = float(metadata["frame_interval_s"])
    active_frames = int(metadata["active_frame_count"])
    rows: list[dict[str, Any]] = []
    for particle_id, trajectory in unit.groupby("global_particle_id", sort=True):
        frames = np.sort(trajectory["frame"].unique())
        first_frame = int(frames[0])
        last_frame = int(frames[-1])
        rows.append(
            {
                "condition": metadata["condition"],
                "source_id": metadata["source_id"],
                "local_cell_id": metadata["local_cell_id"],
                "global_cell_id": metadata["global_cell_id"],
                "channel": int(metadata["channel"]),
                "original_particle_id": trajectory["original_particle_id"].iloc[0],
                "global_particle_id": particle_id,
                "n_trajectories_for_particle": 1,
                "n_detections": int(len(trajectory)),
                "n_unique_frames": int(len(frames)),
                "first_frame": first_frame,
                "last_frame": last_frame,
                "span_frames": last_frame - first_frame + 1,
                "n_missing_frames": last_frame - first_frame + 1 - len(frames),
                "has_frame_gaps": bool(
                    len(frames) > 1 and np.any(np.diff(frames) != 1)
                ),
                "elapsed_duration_s": (last_frame - first_frame) * dt,
                "movie_fraction_observed": len(frames) / active_frames,
                "maximum_supported_lag_frames": last_frame - first_frame,
                "maximum_supported_lag_s": (last_frame - first_frame) * dt,
            }
        )
    return pd.DataFrame(rows)


def _trajectory_filter_summary(
    tracking: pd.DataFrame,
    *,
    minimum_trajectory_frames: int,
    reject_trajectory_gaps: bool,
) -> pd.DataFrame:
    """Describe and classify every trajectory before MSD calculations."""
    rows: list[dict[str, Any]] = []
    for particle_id, trajectory in tracking.groupby("global_particle_id", sort=True):
        metadata = trajectory.iloc[0]
        frames = np.sort(trajectory["frame"].dropna().astype(int).unique())
        first_frame = int(frames[0])
        last_frame = int(frames[-1])
        span_frames = last_frame - first_frame + 1
        n_unique_frames = int(len(frames))
        n_missing_frames = int(span_frames - n_unique_frames)
        has_frame_gaps = bool(
            len(frames) > 1 and np.any(np.diff(frames) != 1)
        )
        too_short = n_unique_frames < int(minimum_trajectory_frames)
        rejected_for_gaps = bool(reject_trajectory_gaps and has_frame_gaps)
        reasons = []
        if too_short:
            reasons.append(f"fewer_than_{int(minimum_trajectory_frames)}_frames")
        if rejected_for_gaps:
            reasons.append("frame_gaps")
        rows.append(
            {
                "condition": metadata["condition"],
                "source_id": metadata["source_id"],
                "local_cell_id": metadata["local_cell_id"],
                "global_cell_id": metadata["global_cell_id"],
                "channel": int(metadata["channel"]),
                "original_particle_id": metadata["original_particle_id"],
                "global_particle_id": particle_id,
                "n_detections": int(len(trajectory)),
                "n_unique_frames": n_unique_frames,
                "first_frame": first_frame,
                "last_frame": last_frame,
                "span_frames": span_frames,
                "n_missing_frames": n_missing_frames,
                "has_frame_gaps": has_frame_gaps,
                "eligible": not (too_short or rejected_for_gaps),
                "exclusion_reason": ";".join(reasons),
            }
        )
    return pd.DataFrame(rows)


def _individual_msd_for_unit(
    unit: pd.DataFrame,
    *,
    max_lag_frames: int,
    prepared_positions: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Calculate particle-level MSD curves and lag summaries for one unit."""
    metadata = unit.iloc[0]
    dt = float(metadata["frame_interval_s"])
    mpp = float(metadata["xy_um_per_px"])
    is_3d = bool(metadata["is_3d"])
    position_columns = ["x", "y"]
    if prepared_positions is None:
        columns = ["frame", "x", "y", "global_particle_id"]
        if is_3d:
            columns.append("z")
            position_columns.append("z")
        prepared = unit[columns].copy()
        prepared = prepared.rename(columns={"global_particle_id": "particle"})
        if is_3d:
            prepared["z"] *= float(metadata["z_um_per_px"]) / mpp
    else:
        prepared = prepared_positions.copy()
        if is_3d:
            position_columns.append("z")

    individual = tp.imsd(
        prepared,
        mpp=mpp,
        fps=1.0 / dt,
        max_lagtime=max_lag_frames,
        pos_columns=position_columns,
    )
    if isinstance(individual, pd.Series):
        individual = individual.to_frame()
    individual.index = pd.Index(
        individual.index.to_numpy(dtype=float), name="lag_seconds"
    )

    n = individual.notna().sum(axis=1).astype(int)
    mean = individual.mean(axis=1, skipna=True)
    sd = individual.std(axis=1, skipna=True, ddof=1)
    summary = pd.DataFrame(
        {
            "lag_seconds": individual.index.to_numpy(dtype=float),
            "lag_frames": np.rint(individual.index.to_numpy(dtype=float) / dt).astype(
                int
            ),
            "particle_mean_msd_um2": mean.to_numpy(dtype=float),
            "particle_sd_msd_um2": sd.to_numpy(dtype=float),
            "particle_sem_msd_um2": (sd / np.sqrt(n.where(n > 0, np.nan))).to_numpy(
                dtype=float
            ),
            "n_contributing_trajectories": n.to_numpy(dtype=int),
        }
    )

    long = (
        individual.rename_axis(index="lag_seconds", columns="global_particle_id")
        .reset_index()
        .melt(
            id_vars="lag_seconds",
            var_name="global_particle_id",
            value_name="msd_um2",
        )
        .dropna(subset=["msd_um2"])
        .reset_index(drop=True)
    )
    long["lag_frames"] = np.rint(long["lag_seconds"] / dt).astype(int)
    for column in (
        "condition",
        "source_id",
        "local_cell_id",
        "global_cell_id",
        "channel",
    ):
        long[column] = metadata[column]
        summary[column] = metadata[column]
    return long, summary


def _fit_result_row(
    result: Any | None,
    *,
    model: str,
    metadata: Mapping[str, Any],
    n_trajectories: int,
    requested_fit_lag_s: float,
    requested_fit_points: int,
    failure_reason: str = "",
) -> dict[str, Any]:
    """Convert a MicroLive fit result into one tidy summary row."""
    mpp = float(metadata["xy_um_per_px"])
    base = {
        "condition": metadata["condition"],
        "source_id": metadata["source_id"],
        "local_cell_id": metadata["local_cell_id"],
        "global_cell_id": metadata["global_cell_id"],
        "channel": int(metadata["channel"]),
        "model": model,
        "dimensions": 3 if bool(metadata["is_3d"]) else 2,
        "n_trajectories": int(n_trajectories),
        "requested_fit_lag_s": float(requested_fit_lag_s),
        "fit_points_requested": int(requested_fit_points),
        "fit_points_used": 0,
        "fit_lag_start_s": np.nan,
        "fit_lag_end_s": np.nan,
        "D_um2_s": np.nan,
        "D_se_um2_s": np.nan,
        "D_ci95_low_um2_s": np.nan,
        "D_ci95_high_um2_s": np.nan,
        "D_px2_s": np.nan,
        "K_alpha_um2_s_alpha": np.nan,
        "K_alpha_se_um2_s_alpha": np.nan,
        "K_alpha_ci95_low_um2_s_alpha": np.nan,
        "K_alpha_ci95_high_um2_s_alpha": np.nan,
        "K_alpha_px2_s_alpha": np.nan,
        "alpha": np.nan,
        "alpha_se": np.nan,
        "alpha_ci95_low": np.nan,
        "alpha_ci95_high": np.nan,
        "offset_um2": np.nan,
        "r_squared": np.nan,
        "pseudo_r_squared": np.nan,
        "rss": np.nan,
        "rmse": np.nan,
        "aicc": np.nan,
        "classification": "Indeterminate",
        "optimizer_success": False,
        "fit_valid": False,
        "at_parameter_bound": False,
        "weakly_identified": False,
        "fit_status": failure_reason,
        "xy_um_per_px": mpp,
        "z_um_per_px": metadata["z_um_per_px"],
        "frame_interval_s": metadata["frame_interval_s"],
        "unit_D": "µm²/s",
        "unit_K_alpha": "µm²/s^alpha",
    }
    if result is None:
        return base

    values = result.parameter_values
    errors = result.parameter_standard_errors
    ci_d = result.parameter_ci95.get("D", (np.nan, np.nan))
    ci_k = result.parameter_ci95.get("K_alpha", (np.nan, np.nan))
    ci_alpha = result.parameter_ci95.get("alpha", (np.nan, np.nan))
    status_parts = [result.status_message] if result.status_message else []
    status_parts.extend(result.warnings)
    base.update(
        {
            "fit_points_used": result.n_fit_points,
            "fit_lag_start_s": result.fit_lag_start_s,
            "fit_lag_end_s": result.fit_lag_end_s,
            "D_um2_s": values.get("D", np.nan),
            "D_se_um2_s": errors.get("D", np.nan),
            "D_ci95_low_um2_s": ci_d[0],
            "D_ci95_high_um2_s": ci_d[1],
            "D_px2_s": (values["D"] / (mpp**2) if "D" in values else np.nan),
            "K_alpha_um2_s_alpha": values.get("K_alpha", np.nan),
            "K_alpha_se_um2_s_alpha": errors.get("K_alpha", np.nan),
            "K_alpha_ci95_low_um2_s_alpha": ci_k[0],
            "K_alpha_ci95_high_um2_s_alpha": ci_k[1],
            "K_alpha_px2_s_alpha": (
                values["K_alpha"] / (mpp**2) if "K_alpha" in values else np.nan
            ),
            "alpha": values.get("alpha", np.nan),
            "alpha_se": errors.get("alpha", np.nan),
            "alpha_ci95_low": ci_alpha[0],
            "alpha_ci95_high": ci_alpha[1],
            "offset_um2": result.offset_um2,
            "r_squared": result.r_squared,
            "pseudo_r_squared": result.pseudo_r_squared,
            "rss": result.rss,
            "rmse": result.rmse,
            "aicc": result.aicc,
            "classification": result.classification,
            "optimizer_success": bool(result.success),
            "fit_valid": bool(result.fit_reliable),
            "at_parameter_bound": bool(result.at_parameter_bound),
            "weakly_identified": bool(result.weakly_identified),
            "fit_status": "; ".join(dict.fromkeys(filter(None, status_parts))),
        }
    )
    return base


def _empty_unit_summary(
    metadata: Mapping[str, Any],
    n_trajectories: int,
    reason: str,
) -> dict[str, Any]:
    """Create an excluded cell summary row for a failed analysis unit."""
    return {
        "condition": metadata["condition"],
        "source_id": metadata["source_id"],
        "local_cell_id": metadata["local_cell_id"],
        "global_cell_id": metadata["global_cell_id"],
        "channel": int(metadata["channel"]),
        "n_trajectories": int(n_trajectories),
        "median_trajectory_frames": np.nan,
        "D_um2_s": np.nan,
        "D_fit_valid": False,
        "alpha": np.nan,
        "K_alpha_um2_s_alpha": np.nan,
        "anomalous_fit_valid": False,
        "classification": "Indeterminate",
        "exclusion_reason": reason,
    }


def _analyze_tracking_results(
    results_root: str | Path,
    *,
    maximum_fit_lag_seconds: float = 25.0,
    maximum_computed_lag_frames: int | None = None,
    minimum_trajectory_frames: int = 1,
    reject_trajectory_gaps: bool = False,
    channels: Sequence[int] | None = None,
    minimum_contributors_per_lag: int = 2,
    remove_drift: bool = False,
    strict: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Analyze every source cell and tracking channel.

    Normal and anomalous fits come directly from ``mi.ParticleMotion``.
    ``maximum_fit_lag_seconds`` controls the shared early-lag fit interval for
    both models. ``minimum_trajectory_frames`` rejects trajectories with fewer
    unique observed frames before MSD and spot-metric calculations. When
    ``reject_trajectory_gaps`` is true, a trajectory is also rejected if any
    frame is missing between its first and last observation. The returned
    dictionary contains tidy tables plus the MicroLive result objects needed
    to render fits without recalculation.
    """
    if not np.isfinite(maximum_fit_lag_seconds) or maximum_fit_lag_seconds <= 0:
        raise ValueError("maximum_fit_lag_seconds must be positive.")
    if maximum_computed_lag_frames is not None and maximum_computed_lag_frames < 1:
        raise ValueError("maximum_computed_lag_frames must be at least one.")
    if (
        isinstance(minimum_trajectory_frames, (bool, np.bool_))
        or int(minimum_trajectory_frames) != minimum_trajectory_frames
        or int(minimum_trajectory_frames) < 1
    ):
        raise ValueError("minimum_trajectory_frames must be a positive integer.")
    if not isinstance(reject_trajectory_gaps, (bool, np.bool_)):
        raise ValueError("reject_trajectory_gaps must be a boolean.")
    if minimum_contributors_per_lag < 1:
        raise ValueError("minimum_contributors_per_lag must be at least one.")

    inventory, tracking = _load_results(
        results_root,
        strict=strict,
        verbose=verbose,
    )
    detected_channels = tuple(sorted(tracking["channel"].unique().tolist()))
    selected_channels = (
        detected_channels
        if channels is None
        else tuple(dict.fromkeys(int(channel) for channel in channels))
    )
    if not selected_channels:
        raise ValueError("At least one channel must be selected.")
    absent_globally = sorted(set(selected_channels).difference(detected_channels))
    if absent_globally:
        raise ValueError(f"Requested channels are absent: {absent_globally}.")

    trajectory_filter_summary = _trajectory_filter_summary(
        tracking,
        minimum_trajectory_frames=int(minimum_trajectory_frames),
        reject_trajectory_gaps=bool(reject_trajectory_gaps),
    )
    retained_particle_ids = set(
        trajectory_filter_summary.loc[
            trajectory_filter_summary["eligible"], "global_particle_id"
        ]
    )
    filtered_tracking = tracking.loc[
        tracking["global_particle_id"].isin(retained_particle_ids)
    ].copy()
    rejected_short_count = int(
        (
            trajectory_filter_summary["n_unique_frames"]
            < int(minimum_trajectory_frames)
        ).sum()
    )
    rejected_gap_count = int(
        (
            (
                trajectory_filter_summary["n_unique_frames"]
                >= int(minimum_trajectory_frames)
            )
            & trajectory_filter_summary["has_frame_gaps"]
            & bool(reject_trajectory_gaps)
        ).sum()
    )
    if filtered_tracking.empty:
        raise ValueError(
            "No trajectories remain after applying trajectory acceptance: "
            f"minimum_trajectory_frames={minimum_trajectory_frames}, "
            f"reject_trajectory_gaps={bool(reject_trajectory_gaps)}."
        )

    trajectory_parts: list[pd.DataFrame] = []
    particle_msd_parts: list[pd.DataFrame] = []
    cell_msd_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    cell_summary_rows: list[dict[str, Any]] = []
    unit_results: list[dict[str, Any]] = []

    cell_groups = tracking.groupby(
        ["condition", "source_id", "global_cell_id"], sort=True
    )
    for (_, _, global_cell_id), cell in cell_groups:
        for channel in selected_channels:
            unit = cell[cell["channel"] == channel].copy()
            if unit.empty:
                template = cell.iloc[0].to_dict()
                template["channel"] = channel
                template["local_cell_id"] = cell["local_cell_id"].iloc[0]
                reason = f"Channel {channel} is absent from this cell."
                fit_points = int(
                    math.floor(
                        maximum_fit_lag_seconds / float(template["frame_interval_s"])
                        + 1e-9
                    )
                )
                for model in ("normal", "anomalous"):
                    fit_rows.append(
                        _fit_result_row(
                            None,
                            model=model,
                            metadata=template,
                            n_trajectories=0,
                            requested_fit_lag_s=maximum_fit_lag_seconds,
                            requested_fit_points=fit_points,
                            failure_reason=reason,
                        )
                    )
                cell_summary_rows.append(_empty_unit_summary(template, 0, reason))
                continue

            unit = unit.loc[
                unit["global_particle_id"].isin(retained_particle_ids)
            ].copy()
            if unit.empty:
                template = cell.iloc[0].to_dict()
                template["channel"] = channel
                template["local_cell_id"] = cell["local_cell_id"].iloc[0]
                reason = (
                    f"All trajectories in channel {channel} failed acceptance "
                    f"(minimum {minimum_trajectory_frames} unique frames; "
                    f"reject gaps={bool(reject_trajectory_gaps)})."
                )
                fit_points = int(
                    math.floor(
                        maximum_fit_lag_seconds / float(template["frame_interval_s"])
                        + 1e-9
                    )
                )
                for model in ("normal", "anomalous"):
                    fit_rows.append(
                        _fit_result_row(
                            None,
                            model=model,
                            metadata=template,
                            n_trajectories=0,
                            requested_fit_lag_s=maximum_fit_lag_seconds,
                            requested_fit_points=fit_points,
                            failure_reason=reason,
                        )
                    )
                cell_summary_rows.append(_empty_unit_summary(template, 0, reason))
                if verbose:
                    print(f"Warning: {global_cell_id}, channel {channel}: {reason}")
                continue

            metadata = unit.iloc[0].to_dict()
            trajectory_summary = _trajectory_summary_for_unit(unit)
            trajectory_parts.append(trajectory_summary)
            n_trajectories = int(unit["global_particle_id"].nunique())
            median_track_frames = float(trajectory_summary["n_unique_frames"].median())
            dt = float(metadata["frame_interval_s"])
            frame_span = int(unit["frame"].max() - unit["frame"].min())
            computed_lag_frames = (
                frame_span
                if maximum_computed_lag_frames is None
                else min(int(maximum_computed_lag_frames), frame_span)
            )
            requested_fit_points = int(math.floor(maximum_fit_lag_seconds / dt + 1e-9))
            reason = ""
            try:
                if computed_lag_frames < 1:
                    raise ValueError(
                        "The cell/channel contains fewer than two movie frames."
                    )
                if requested_fit_points < 2:
                    raise ValueError(
                        f"The requested fit lag ({maximum_fit_lag_seconds:g} s) "
                        f"contains fewer than two positive lags at Δt={dt:g} s."
                    )

                motion_input = unit.copy()
                motion_input["unique_particle"] = motion_input["global_particle_id"]
                motion = mi.ParticleMotion(
                    trackpy_dataframe=motion_input,
                    microns_per_pixel=float(metadata["xy_um_per_px"]),
                    step_size_in_sec=dt,
                    max_lagtime=computed_lag_frames,
                    show_plot=False,
                    remove_drift=remove_drift,
                    spot_type=int(channel),
                    max_fit_points=requested_fit_points,
                    is_3d=bool(metadata["is_3d"]),
                    microns_per_pixel_z=(
                        float(metadata["z_um_per_px"])
                        if bool(metadata["is_3d"])
                        else None
                    ),
                )
                (
                    _,
                    _,
                    ensemble_msd,
                    _,
                    _,
                    _,
                    msd_positions,
                ) = motion.calculate_msd()

                # Use the same coordinates for the particle and ensemble MSDs.
                # ``msd_positions`` is drift-corrected when requested.
                particle_long, particle_summary = _individual_msd_for_unit(
                    unit,
                    max_lag_frames=computed_lag_frames,
                    prepared_positions=msd_positions,
                )
                particle_long = particle_long.merge(
                    trajectory_summary[
                        [
                            "global_particle_id",
                            "n_detections",
                            "n_unique_frames",
                            "first_frame",
                            "last_frame",
                            "elapsed_duration_s",
                        ]
                    ],
                    on="global_particle_id",
                    how="left",
                    validate="many_to_one",
                )

                ensemble = pd.DataFrame(
                    {
                        "lag_seconds": ensemble_msd.index.to_numpy(dtype=float),
                        "lag_frames": np.rint(
                            ensemble_msd.index.to_numpy(dtype=float) / dt
                        ).astype(int),
                        "pair_weighted_ensemble_msd_um2": ensemble_msd.to_numpy(
                            dtype=float
                        ),
                    }
                )
                for column in (
                    "condition",
                    "source_id",
                    "local_cell_id",
                    "global_cell_id",
                    "channel",
                ):
                    ensemble[column] = metadata[column]
                cell_curve = ensemble.merge(
                    particle_summary.drop(columns="lag_seconds"),
                    on=[
                        "condition",
                        "source_id",
                        "local_cell_id",
                        "global_cell_id",
                        "channel",
                        "lag_frames",
                    ],
                    how="outer",
                    validate="one_to_one",
                )
                # Trackpy can represent the same physical lag as, for example,
                # 10.0 or 9.999999999999998 in different estimators. Join on
                # integer frame lag, then reconstruct the exact physical time.
                cell_curve["lag_seconds"] = (
                    cell_curve["lag_frames"].to_numpy(dtype=float) * dt
                )
                cell_curve = cell_curve.sort_values("lag_frames")
                cell_curve["minimum_contributors_met"] = (
                    cell_curve["n_contributing_trajectories"].fillna(0)
                    >= minimum_contributors_per_lag
                )

                fit_data = cell_curve.loc[
                    cell_curve["minimum_contributors_met"]
                    & (cell_curve["lag_seconds"] > 0)
                    & (cell_curve["lag_seconds"] <= maximum_fit_lag_seconds),
                    ["lag_seconds", "pair_weighted_ensemble_msd_um2"],
                ].dropna()
                if len(fit_data) < 2:
                    raise ValueError(
                        "Fewer than two MSD lag points meet the fit interval and "
                        f"minimum contributor threshold ({minimum_contributors_per_lag})."
                    )
                normal_fit, anomalous_fit = motion.fit_msd_models(
                    fit_data["lag_seconds"].to_numpy(dtype=float),
                    fit_data["pair_weighted_ensemble_msd_um2"].to_numpy(dtype=float),
                )
                fit_times = normal_fit.fit_times_s
                cell_curve["is_fit_lag"] = cell_curve["lag_seconds"].apply(
                    lambda value: bool(
                        np.any(
                            np.isclose(
                                float(value),
                                np.asarray(fit_times, dtype=float),
                                rtol=1e-7,
                                atol=max(1e-10, dt * 1e-8),
                            )
                        )
                    )
                )
                particle_msd_parts.append(particle_long)
                cell_msd_parts.append(cell_curve)
                fit_rows.extend(
                    [
                        _fit_result_row(
                            normal_fit,
                            model="normal",
                            metadata=metadata,
                            n_trajectories=n_trajectories,
                            requested_fit_lag_s=maximum_fit_lag_seconds,
                            requested_fit_points=requested_fit_points,
                        ),
                        _fit_result_row(
                            anomalous_fit,
                            model="anomalous",
                            metadata=metadata,
                            n_trajectories=n_trajectories,
                            requested_fit_lag_s=maximum_fit_lag_seconds,
                            requested_fit_points=requested_fit_points,
                        ),
                    ]
                )
                cell_summary_rows.append(
                    {
                        "condition": metadata["condition"],
                        "source_id": metadata["source_id"],
                        "local_cell_id": metadata["local_cell_id"],
                        "global_cell_id": global_cell_id,
                        "channel": int(channel),
                        "n_trajectories": n_trajectories,
                        "median_trajectory_frames": median_track_frames,
                        "D_um2_s": normal_fit.parameter_values.get("D", np.nan),
                        "D_fit_valid": bool(normal_fit.fit_reliable),
                        "alpha": anomalous_fit.parameter_values.get("alpha", np.nan),
                        "K_alpha_um2_s_alpha": anomalous_fit.parameter_values.get(
                            "K_alpha", np.nan
                        ),
                        "anomalous_fit_valid": bool(anomalous_fit.fit_reliable),
                        "classification": anomalous_fit.classification,
                        "exclusion_reason": "",
                    }
                )
                unit_results.append(
                    {
                        "condition": metadata["condition"],
                        "source_id": metadata["source_id"],
                        "local_cell_id": metadata["local_cell_id"],
                        "global_cell_id": global_cell_id,
                        "channel": int(channel),
                        "n_trajectories": n_trajectories,
                        "normal_fit": normal_fit,
                        "anomalous_fit": anomalous_fit,
                    }
                )
            except Exception as exc:
                reason = str(exc)
                if strict:
                    raise RuntimeError(
                        f"MSD analysis failed for {global_cell_id}, channel {channel}: "
                        f"{reason}"
                    ) from exc
                for model in ("normal", "anomalous"):
                    fit_rows.append(
                        _fit_result_row(
                            None,
                            model=model,
                            metadata=metadata,
                            n_trajectories=n_trajectories,
                            requested_fit_lag_s=maximum_fit_lag_seconds,
                            requested_fit_points=requested_fit_points,
                            failure_reason=reason,
                        )
                    )
                failed = _empty_unit_summary(metadata, n_trajectories, reason)
                failed["median_trajectory_frames"] = median_track_frames
                cell_summary_rows.append(failed)
                if verbose:
                    print(f"Warning: {global_cell_id}, channel {channel}: {reason}")

    trajectory_summary = (
        pd.concat(trajectory_parts, ignore_index=True, sort=False)
        if trajectory_parts
        else pd.DataFrame()
    )
    particle_msd = (
        pd.concat(particle_msd_parts, ignore_index=True, sort=False)
        if particle_msd_parts
        else pd.DataFrame()
    )
    cell_msd = (
        pd.concat(cell_msd_parts, ignore_index=True, sort=False)
        if cell_msd_parts
        else pd.DataFrame()
    )
    fit_summary = pd.DataFrame(fit_rows)
    cell_summary = pd.DataFrame(cell_summary_rows)

    configuration = {
        "results_root": str(Path(results_root).expanduser().resolve()),
        "maximum_fit_lag_seconds": float(maximum_fit_lag_seconds),
        "maximum_computed_lag_frames": maximum_computed_lag_frames,
        "minimum_trajectory_frames": int(minimum_trajectory_frames),
        "reject_trajectory_gaps": bool(reject_trajectory_gaps),
        "channels": selected_channels,
        "minimum_contributors_per_lag": int(minimum_contributors_per_lag),
        "remove_drift": bool(remove_drift),
        "strict_validation": bool(strict),
    }
    analysis = {
        "configuration": configuration,
        "dataset_inventory": inventory,
        "trajectory_filter_summary": trajectory_filter_summary,
        "tracking_data": filtered_tracking,
        "trajectory_summary": trajectory_summary,
        "particle_msd": particle_msd,
        "cell_msd": cell_msd,
        "fit_summary": fit_summary,
        "cell_summary": cell_summary,
        "unit_results": unit_results,
    }

    if verbose:
        counts = (
            trajectory_summary.groupby("channel")["global_particle_id"].nunique()
            if not trajectory_summary.empty
            else pd.Series(dtype=int)
        )
        print(
            f"Analyzed {len(unit_results)} cell/channel units. "
            f"Trajectory totals by channel: {counts.to_dict()}"
        )
        if rejected_short_count:
            print(
                f"Rejected {rejected_short_count} trajectories shorter than "
                f"{minimum_trajectory_frames} frames."
            )
        if rejected_gap_count:
            print(
                f"Rejected {rejected_gap_count} additional trajectories with "
                "internal frame gaps."
            )
    return analysis


def _save_base_tables(
    analysis: Mapping[str, Any],
    output_directory: str | Path,
) -> dict[str, Path]:
    """Save the tidy numerical outputs; original result folders stay read-only."""
    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for key, filename in _BASE_TABLE_FILENAMES.items():
        table = analysis.get(key)
        if not isinstance(table, pd.DataFrame):
            continue
        path = output / filename
        table.to_csv(path, index=False)
        written[key] = path
    return written


def _publication_font() -> str:
    """Select Arial when available and otherwise use Matplotlib's sans font."""
    try:
        font_manager.findfont(
            font_manager.FontProperties(family=["Arial"]),
            fallback_to_default=False,
        )
        return "Arial"
    except ValueError:
        warnings.warn(
            "Arial is not installed or is not visible to Matplotlib; "
            "falling back to DejaVu Sans.",
            RuntimeWarning,
            stacklevel=2,
        )
        return "DejaVu Sans"


def _publication_style(
    overrides: Mapping[str, Any] | None = None,
) -> mpl.rc_context:
    """Return the shared white-background, no-grid publication style context."""
    font_name = _publication_font()
    settings = {
        "font.family": "sans-serif",
        "font.sans-serif": [font_name],
        "font.size": 18,
        "text.color": "black",
        "axes.labelsize": 24,
        "axes.labelweight": "semibold",
        "axes.labelcolor": "black",
        "axes.titlesize": 22,
        "axes.titleweight": "semibold",
        "axes.titlecolor": "black",
        "axes.edgecolor": "black",
        "axes.linewidth": 1.6,
        "axes.facecolor": "white",
        "axes.grid": False,
        "figure.facecolor": "white",
        "figure.edgecolor": "white",
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "savefig.transparent": False,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "xtick.color": "black",
        "ytick.color": "black",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
        "legend.fontsize": 15,
        "legend.frameon": True,
        "legend.edgecolor": "black",
        "legend.facecolor": "white",
        "legend.framealpha": 0.95,
        "svg.fonttype": "none",
    }
    if overrides:
        settings.update(dict(overrides))
    return mpl.rc_context(settings)


# Paired-channel summaries, statistics, and figures


def _normalize_selection(
    available: Iterable[int],
    selected: Sequence[int] | int | None,
) -> list[int]:
    """Validate and deduplicate the requested channel selection."""
    available_list = sorted(int(value) for value in available)
    if selected is None:
        return available_list
    if isinstance(selected, (int, np.integer)):
        values = [int(selected)]
    else:
        values = [int(value) for value in selected]
    values = list(dict.fromkeys(values))
    missing = sorted(set(values).difference(available_list))
    if missing:
        raise ValueError(f"Requested channels are absent: {missing}.")
    if not values:
        raise ValueError("At least one channel must be selected.")
    return values


def _resolve_channel_condition_labels(
    channels: Sequence[int],
    labels: Mapping[int, str] | None,
) -> dict[int, str]:
    """Resolve unique display labels for selected channels."""
    provided = {} if labels is None else {int(k): str(v) for k, v in labels.items()}
    resolved: dict[int, str] = {}
    for channel in channels:
        label = provided.get(int(channel), f"Channel {int(channel)}").strip()
        if not label:
            raise ValueError(f"The label for channel {channel} is empty.")
        resolved[int(channel)] = label
    if len(set(resolved.values())) != len(resolved):
        raise ValueError("Channel-condition labels must be unique.")
    return resolved


def _resolve_measurement_channel_map(
    channels: Sequence[int],
    mapping: Mapping[int, int] | None,
) -> dict[int, int]:
    """Resolve measurement-channel assignments for selected channels."""
    provided = {} if mapping is None else {int(k): int(v) for k, v in mapping.items()}
    return {
        int(channel): int(provided.get(int(channel), int(channel)))
        for channel in channels
    }


def _validate_reducer(name: str) -> str:
    """Validate and normalize a trajectory or cell aggregation reducer."""
    normalized = str(name).strip().lower()
    if normalized not in {"mean", "median"}:
        raise ValueError("Reducer must be 'mean' or 'median'.")
    return normalized


def _channel_encoding(channels: Sequence[int]) -> dict[int, tuple[str, str]]:
    """Assign repeatable colors and markers to selected channels."""
    return {
        int(channel): (
            _CHANNEL_COLORS[index % len(_CHANNEL_COLORS)],
            _CHANNEL_MARKERS[index % len(_CHANNEL_MARKERS)],
        )
        for index, channel in enumerate(channels)
    }


def _add_channel_condition_column(
    table: pd.DataFrame,
    labels: Mapping[int, str],
) -> pd.DataFrame:
    """Add human-readable channel-condition labels to a table copy."""
    copied = table.copy()
    if "channel" in copied.columns:
        copied["channel_condition"] = copied["channel"].map(labels)
    return copied


def _save_figure(
    figure: plt.Figure,
    *,
    output_directory: str | Path,
    stem: str,
    dpi: int,
) -> tuple[Path, Path]:
    """Save a figure as PNG and SVG using a sanitized filename stem."""
    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("_")
    png = output / f"{safe_stem}.png"
    svg = output / f"{safe_stem}.svg"
    figure.savefig(
        png,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="white",
        transparent=False,
    )
    figure.savefig(
        svg,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="white",
        transparent=False,
    )
    return png, svg


def fit_lag_seconds_from_frames(
    results_root: str | Path,
    maximum_fit_lag_frames: int,
) -> float:
    """Convert a user-selected lag-frame count using result metadata.

    Only metadata files are read. All retained result folders must report the
    same positive frame interval so that one frame-based fit setting has an
    unambiguous physical-time meaning.
    """
    if (
        isinstance(maximum_fit_lag_frames, (bool, np.bool_))
        or int(maximum_fit_lag_frames) != maximum_fit_lag_frames
        or int(maximum_fit_lag_frames) < 1
    ):
        raise ValueError("maximum_fit_lag_frames must be a positive integer.")
    root = Path(results_root).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Results directory not found: {root}")
    result_folders = sorted(
        folder
        for folder in root.iterdir()
        if folder.is_dir() and folder.name.startswith("results_")
    )
    if not result_folders:
        raise ValueError(f"No results_* folders found in {root}.")

    intervals: list[float] = []
    for folder in result_folders:
        metadata_files = sorted(
            path
            for path in folder.glob("Metadata_*.txt")
            if not path.name.startswith("._")
        )
        if len(metadata_files) != 1:
            raise ValueError(
                f"{folder.name} must contain exactly one Metadata_*.txt file; "
                f"found {len(metadata_files)}."
            )
        metadata = _read_microlive_metadata(metadata_files[0])
        intervals.append(float(metadata["frame_interval_s"]))

    reference = intervals[0]
    if not np.allclose(intervals, reference, rtol=1e-9, atol=1e-12):
        unique = sorted(set(intervals))
        raise ValueError(
            "A frame-based shared fit interval requires one common frame "
            f"interval; found {unique} s/frame."
        )
    return float(int(maximum_fit_lag_frames) * reference)


def summarize_spot_metrics(
    tracking_data: pd.DataFrame,
    *,
    channels: Sequence[int] | int | None = None,
    channel_condition_labels: Mapping[int, str] | None = None,
    measurement_channel_map: Mapping[int, int] | None = None,
    trajectory_reducer: str = "mean",
    cell_reducer: str = "mean",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate matched-channel spot measurements trajectory -> cell.

    Long trajectories do not receive extra cell-level weight: detections are
    summarized once per trajectory before trajectories are summarized per cell.
    """
    required = {
        "condition",
        "source_id",
        "local_cell_id",
        "global_cell_id",
        "channel",
        "global_particle_id",
        "frame",
        "xy_um_per_px",
        "snr_method",
    }
    missing = sorted(required.difference(tracking_data.columns))
    if missing:
        raise ValueError(f"tracking_data is missing required columns: {missing}.")
    if tracking_data.empty:
        raise ValueError("tracking_data is empty.")

    duplicate = tracking_data.duplicated(["global_particle_id", "frame"], keep=False)
    if duplicate.any():
        raise ValueError(
            "tracking_data contains duplicate (global_particle_id, frame) rows."
        )

    particle_channels = tracking_data.groupby("global_particle_id")["channel"].nunique()
    if (particle_channels > 1).any():
        raise ValueError("At least one global particle maps to multiple channels.")

    trajectory_reducer = _validate_reducer(trajectory_reducer)
    cell_reducer = _validate_reducer(cell_reducer)
    available_channels = sorted(tracking_data["channel"].dropna().astype(int).unique())
    selected_channels = _normalize_selection(available_channels, channels)
    labels = _resolve_channel_condition_labels(
        selected_channels, channel_condition_labels
    )
    measurement_map = _resolve_measurement_channel_map(
        selected_channels, measurement_channel_map
    )

    parts: list[pd.DataFrame] = []
    for channel in selected_channels:
        measurement_channel = measurement_map[channel]
        intensity_column = f"spot_int_ch_{measurement_channel}"
        snr_column = f"snr_ch_{measurement_channel}"
        size_column = f"spot_size_ch_{measurement_channel}"
        absent = [
            column
            for column in (intensity_column, snr_column, size_column)
            if column not in tracking_data.columns
        ]
        if absent:
            raise ValueError(
                f"Tracking channel {channel} is missing matched measurement "
                f"columns: {absent}."
            )

        selected = tracking_data.loc[tracking_data["channel"] == channel].copy()
        snr_methods = set(
            selected["snr_method"]
            .fillna("missing")
            .astype(str)
            .str.strip()
            .str.lower()
            .unique()
        )
        if snr_methods != {"disk_doughnut"}:
            raise ValueError(
                f"Tracking channel {channel} must use MicroLive's "
                f"disk_doughnut SNR method; found {sorted(snr_methods)}."
            )
        calibration = pd.to_numeric(selected["xy_um_per_px"], errors="coerce").to_numpy(
            dtype=float
        )
        if not np.isfinite(calibration).all() or np.any(calibration <= 0):
            raise ValueError(f"Tracking channel {channel} has invalid XY calibration.")

        selected["measurement_channel"] = measurement_channel
        selected["channel_condition"] = labels[channel]
        selected["spot_intensity_au"] = pd.to_numeric(
            selected[intensity_column], errors="coerce"
        )
        selected["spot_snr"] = pd.to_numeric(selected[snr_column], errors="coerce")
        selected["spot_fwhm_px"] = pd.to_numeric(selected[size_column], errors="coerce")
        finite_size = selected["spot_fwhm_px"].dropna().to_numpy(dtype=float)
        if finite_size.size and (
            not np.isfinite(finite_size).all() or np.any(finite_size <= 0)
        ):
            raise ValueError(
                f"Tracking channel {channel} contains nonpositive/nonfinite "
                "spot FWHM values."
            )
        selected["spot_fwhm_um"] = selected["spot_fwhm_px"] * pd.to_numeric(
            selected["xy_um_per_px"], errors="coerce"
        )
        selected["spot_fwhm_nm"] = selected["spot_fwhm_um"] * 1000.0
        parts.append(selected)

    prepared = pd.concat(parts, ignore_index=True, sort=False)
    trajectory_keys = [
        "condition",
        "source_id",
        "local_cell_id",
        "global_cell_id",
        "channel",
        "channel_condition",
        "measurement_channel",
        "global_particle_id",
    ]
    trajectory = (
        prepared.groupby(trajectory_keys, as_index=False, sort=True)
        .agg(
            n_detections=("frame", "size"),
            n_unique_frames=("frame", "nunique"),
            first_frame=("frame", "min"),
            last_frame=("frame", "max"),
            trajectory_spot_intensity_au=(
                "spot_intensity_au",
                trajectory_reducer,
            ),
            trajectory_spot_snr=("spot_snr", trajectory_reducer),
            trajectory_spot_fwhm_px=("spot_fwhm_px", trajectory_reducer),
            trajectory_spot_fwhm_um=("spot_fwhm_um", trajectory_reducer),
            trajectory_spot_fwhm_nm=("spot_fwhm_nm", trajectory_reducer),
        )
        .sort_values(["source_id", "local_cell_id", "channel", "global_particle_id"])
        .reset_index(drop=True)
    )

    cell_keys = [
        "condition",
        "source_id",
        "local_cell_id",
        "global_cell_id",
        "channel",
        "channel_condition",
        "measurement_channel",
    ]
    cell = (
        trajectory.groupby(cell_keys, as_index=False, sort=True)
        .agg(
            n_trajectories_spot_metrics=("global_particle_id", "nunique"),
            n_intensity_trajectories=(
                "trajectory_spot_intensity_au",
                "count",
            ),
            n_snr_trajectories=("trajectory_spot_snr", "count"),
            n_fwhm_trajectories=("trajectory_spot_fwhm_um", "count"),
            cell_spot_intensity_au=(
                "trajectory_spot_intensity_au",
                cell_reducer,
            ),
            cell_spot_snr=("trajectory_spot_snr", cell_reducer),
            cell_spot_fwhm_px=("trajectory_spot_fwhm_px", cell_reducer),
            cell_spot_fwhm_um=("trajectory_spot_fwhm_um", cell_reducer),
            cell_spot_fwhm_nm=("trajectory_spot_fwhm_nm", cell_reducer),
        )
        .sort_values(["source_id", "local_cell_id", "channel"])
        .reset_index(drop=True)
    )
    return trajectory, cell


def _metric_valid_rows(data: pd.DataFrame, metric: str) -> pd.DataFrame:
    """Return finite cell-level values that pass metric-specific fit QC."""
    if metric not in data.columns:
        raise ValueError(f"Metric {metric!r} is absent from the cell table.")
    selected = data.copy()
    if metric == "D_um2_s" and "D_fit_valid" in selected.columns:
        selected = selected[selected["D_fit_valid"].fillna(False)]
    elif metric in {"alpha", "K_alpha_um2_s_alpha"}:
        if "anomalous_fit_valid" in selected.columns:
            selected = selected[selected["anomalous_fit_valid"].fillna(False)]
    numeric = pd.to_numeric(selected[metric], errors="coerce")
    selected = selected[np.isfinite(numeric)].copy()
    selected[metric] = numeric[np.isfinite(numeric)]
    return selected


def build_channel_condition_summary(
    cell_analysis_summary: pd.DataFrame,
    *,
    metrics: Sequence[str] | None = None,
    channels: Sequence[int] | int | None = None,
) -> pd.DataFrame:
    """Build long-format descriptive summaries across cell samples."""
    available_channels = sorted(
        cell_analysis_summary["channel"].dropna().astype(int).unique()
    )
    selected_channels = _normalize_selection(available_channels, channels)
    selected_metrics = (
        [
            "D_um2_s",
            "alpha",
            "K_alpha_um2_s_alpha",
            "cell_spot_intensity_au",
            "cell_spot_snr",
            "cell_spot_fwhm_um",
            "cell_spot_fwhm_nm",
            "n_trajectories",
            "median_trajectory_frames",
        ]
        if metrics is None
        else list(metrics)
    )
    rows: list[dict[str, Any]] = []
    for metric in selected_metrics:
        if metric not in cell_analysis_summary.columns:
            continue
        valid = _metric_valid_rows(cell_analysis_summary, metric)
        valid = valid[valid["channel"].isin(selected_channels)]
        for channel in selected_channels:
            group = valid.loc[valid["channel"] == channel]
            values = group[metric].to_numpy(dtype=float)
            n = int(values.size)
            mean = float(np.mean(values)) if n else np.nan
            sd = float(np.std(values, ddof=1)) if n > 1 else np.nan
            sem = float(sd / math.sqrt(n)) if n > 1 else np.nan
            rows.append(
                {
                    "metric": metric,
                    "unit": METRIC_SPECS.get(metric, {}).get("unit", ""),
                    "channel": int(channel),
                    "channel_condition": (
                        group["channel_condition"].iloc[0]
                        if not group.empty
                        else cell_analysis_summary.loc[
                            cell_analysis_summary["channel"] == channel,
                            "channel_condition",
                        ].iloc[0]
                    ),
                    "n_cells": n,
                    "mean": mean,
                    "sd": sd,
                    "sem": sem,
                    "median": float(np.median(values)) if n else np.nan,
                    "q1": float(np.quantile(values, 0.25)) if n else np.nan,
                    "q3": float(np.quantile(values, 0.75)) if n else np.nan,
                    "minimum": float(np.min(values)) if n else np.nan,
                    "maximum": float(np.max(values)) if n else np.nan,
                }
            )
    return pd.DataFrame(rows)


def analyze_final_results(
    results_root: str | Path,
    *,
    maximum_fit_lag_seconds: float = 100.0,
    maximum_computed_lag_frames: int | None = None,
    minimum_trajectory_frames: int = 1,
    reject_trajectory_gaps: bool = False,
    channels: Sequence[int] | None = None,
    minimum_contributors_per_lag: int = 2,
    remove_drift: bool = False,
    channel_condition_labels: Mapping[int, str] | None = None,
    measurement_channel_map: Mapping[int, int] | None = None,
    trajectory_reducer: str = "mean",
    cell_reducer: str = "mean",
    strict: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run MSD analysis with trajectory filtering and add paired summaries."""
    analysis = _analyze_tracking_results(
        results_root,
        maximum_fit_lag_seconds=maximum_fit_lag_seconds,
        maximum_computed_lag_frames=maximum_computed_lag_frames,
        minimum_trajectory_frames=minimum_trajectory_frames,
        reject_trajectory_gaps=reject_trajectory_gaps,
        channels=channels,
        minimum_contributors_per_lag=minimum_contributors_per_lag,
        remove_drift=remove_drift,
        strict=strict,
        verbose=verbose,
    )
    selected_channels = [
        int(channel) for channel in analysis["configuration"]["channels"]
    ]
    labels = _resolve_channel_condition_labels(
        selected_channels, channel_condition_labels
    )
    measurement_map = _resolve_measurement_channel_map(
        selected_channels, measurement_channel_map
    )

    # Source files are pooled by downstream channel summaries. Preserve each
    # source condition in the detailed tables for provenance.
    updated: dict[str, Any] = dict(analysis)
    for key, value in analysis.items():
        if isinstance(value, pd.DataFrame):
            updated[key] = _add_channel_condition_column(value, labels)
    updated_units: list[dict[str, Any]] = []
    for unit in analysis["unit_results"]:
        copied = dict(unit)
        copied["channel_condition"] = labels[int(copied["channel"])]
        updated_units.append(copied)
    updated["unit_results"] = updated_units

    trajectory_spot, cell_spot = summarize_spot_metrics(
        updated["tracking_data"],
        channels=selected_channels,
        channel_condition_labels=labels,
        measurement_channel_map=measurement_map,
        trajectory_reducer=trajectory_reducer,
        cell_reducer=cell_reducer,
    )
    cell_base = updated["cell_summary"].copy()
    merge_keys = [
        "condition",
        "source_id",
        "local_cell_id",
        "global_cell_id",
        "channel",
        "channel_condition",
    ]
    fit_qc_columns = [
        "optimizer_success",
        "at_parameter_bound",
        "weakly_identified",
        "fit_status",
    ]
    fit_qc_parts: list[pd.DataFrame] = []
    for model in ("normal", "anomalous"):
        model_qc = (
            updated["fit_summary"]
            .loc[
                updated["fit_summary"]["model"] == model,
                merge_keys + fit_qc_columns,
            ]
            .copy()
        )
        if model_qc.duplicated(merge_keys).any():
            raise ValueError(
                f"The {model} fit-QC table is not unique per cell/channel."
            )
        model_qc = model_qc.rename(
            columns={column: f"{model}_{column}" for column in fit_qc_columns}
        )
        fit_qc_parts.append(model_qc)
    fit_qc = fit_qc_parts[0].merge(
        fit_qc_parts[1],
        on=merge_keys,
        how="outer",
        validate="one_to_one",
    )
    cell_analysis = cell_base.merge(
        cell_spot,
        on=merge_keys,
        how="left",
        validate="one_to_one",
    )
    cell_analysis = cell_analysis.merge(
        fit_qc,
        on=merge_keys,
        how="left",
        validate="one_to_one",
    )
    cell_analysis["normal_fit_exclusion_reason"] = np.where(
        cell_analysis["D_fit_valid"].fillna(False),
        "",
        cell_analysis["normal_fit_status"].fillna(cell_analysis["exclusion_reason"]),
    )
    cell_analysis["anomalous_fit_exclusion_reason"] = np.where(
        cell_analysis["anomalous_fit_valid"].fillna(False),
        "",
        cell_analysis["anomalous_fit_status"].fillna(cell_analysis["exclusion_reason"]),
    )
    if cell_analysis.duplicated(["global_cell_id", "channel"]).any():
        raise ValueError(
            "The final cell table is not unique on (global_cell_id, channel)."
        )

    configuration = dict(updated["configuration"])
    configuration.update(
        {
            "source_files_pooled": True,
            "statistical_unit": "cell",
            "experimental_condition": "tracking_channel",
            "channel_condition_labels": labels,
            "measurement_channel_map": measurement_map,
            "trajectory_spot_reducer": _validate_reducer(trajectory_reducer),
            "cell_spot_reducer": _validate_reducer(cell_reducer),
            "strict_validation": bool(strict),
        }
    )
    updated["configuration"] = configuration
    updated["trajectory_spot_summary"] = trajectory_spot
    updated["cell_spot_summary"] = cell_spot
    updated["cell_analysis_summary"] = cell_analysis
    updated["channel_condition_summary"] = build_channel_condition_summary(
        cell_analysis,
        channels=selected_channels,
    )
    updated["msd_d_consistency"] = build_channel_msd_summary(
        updated,
        channels=selected_channels,
        center="mean",
        uncertainty="sem",
        paired_cells_only=True,
        require_valid_normal_fit=True,
        normal_fit_max_lag_seconds=maximum_fit_lag_seconds,
    )["D_msd_consistency"]

    if verbose:
        complete_pairs = (
            cell_analysis.groupby("global_cell_id")["channel"]
            .nunique()
            .eq(len(selected_channels))
            .sum()
        )
        print(
            f"Final paired-channel table: {cell_analysis['global_cell_id'].nunique()} "
            f"cells, {len(cell_analysis)} cell/channel rows, "
            f"{int(complete_pairs)} complete channel sets."
        )
    return updated


def save_final_tables(
    analysis: Mapping[str, Any],
    output_directory: str | Path,
) -> dict[str, Path]:
    """Save existing core tables plus final hierarchical tables."""
    written = _save_base_tables(analysis, output_directory)
    output = Path(output_directory).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    for key, filename in _FINAL_TABLE_FILENAMES.items():
        table = analysis.get(key)
        if not isinstance(table, pd.DataFrame) or table.empty:
            continue
        path = output / filename
        table.to_csv(path, index=False)
        written[key] = path
    return written


def _complete_cell_ids(
    cell_summary: pd.DataFrame,
    channels: Sequence[int],
    *,
    require_valid_normal_fit: bool,
) -> list[str]:
    """Return cell identifiers complete across selected channels and QC."""
    selected = cell_summary[cell_summary["channel"].isin(channels)].copy()
    if require_valid_normal_fit and "D_fit_valid" in selected.columns:
        selected = selected[selected["D_fit_valid"].fillna(False)]
    counts = selected.groupby("global_cell_id")["channel"].nunique()
    return sorted(counts[counts == len(channels)].index.tolist())


def _validate_summary_options(center: str, uncertainty: str) -> tuple[str, str]:
    """Validate a statistically coherent center and uncertainty pairing."""
    center = str(center).lower()
    uncertainty = str(uncertainty).lower()
    allowed = {
        "mean": {"sem", "sd", "none"},
        "median": {"iqr", "none"},
    }
    if center not in allowed:
        raise ValueError("center must be 'mean' or 'median'.")
    if uncertainty not in {"sem", "sd", "iqr", "none"}:
        raise ValueError("uncertainty must be 'sem', 'sd', 'iqr', or 'none'.")
    if uncertainty not in allowed[center]:
        choices = ", ".join(sorted(allowed[center]))
        raise ValueError(
            f"uncertainty={uncertainty!r} is not appropriate for center={center!r}; "
            f"choose one of: {choices}."
        )
    return center, uncertainty


def _summary_legend_label(center: str, uncertainty: str) -> str:
    """Format a concise center/uncertainty description for legends."""
    if uncertainty == "none":
        return center
    if uncertainty == "iqr":
        return f"{center} with IQR"
    return f"{center} ± {uncertainty.upper()}"


def _curve_spread(
    values: np.ndarray,
    *,
    center: str,
    uncertainty: str,
) -> dict[str, float]:
    """Summarize the center and uncertainty of one MSD lag distribution."""
    n = int(values.size)
    mean = float(np.mean(values))
    median = float(np.median(values))
    center_value = mean if center == "mean" else median
    sd = float(np.std(values, ddof=1)) if n > 1 else np.nan
    sem = float(sd / math.sqrt(n)) if n > 1 else np.nan
    q1 = float(np.quantile(values, 0.25))
    q3 = float(np.quantile(values, 0.75))
    if uncertainty == "sem":
        lower, upper = center_value - sem, center_value + sem
    elif uncertainty == "sd":
        lower, upper = center_value - sd, center_value + sd
    elif uncertainty == "iqr":
        lower, upper = q1, q3
    else:
        lower = upper = center_value
    return {
        "n_cells": n,
        "mean_msd_um2": mean,
        "median_msd_um2": median,
        "center_msd_um2": center_value,
        "sd_msd_um2": sd,
        "sem_msd_um2": sem,
        "q1_msd_um2": q1,
        "q3_msd_um2": q3,
        "lower_msd_um2": lower,
        "upper_msd_um2": upper,
    }


def build_channel_msd_summary(
    analysis: Mapping[str, Any],
    *,
    channels: Sequence[int] | int | None = None,
    msd_column: str = "pair_weighted_ensemble_msd_um2",
    center: str = "mean",
    uncertainty: str = "sem",
    paired_cells_only: bool = True,
    require_valid_normal_fit: bool = True,
    maximum_display_lag_seconds: float | None = None,
    normal_fit_max_lag_seconds: float | None = None,
) -> dict[str, Any]:
    """Aggregate one cell MSD curve per paired channel condition."""
    center, uncertainty = _validate_summary_options(center, uncertainty)
    if maximum_display_lag_seconds is not None and maximum_display_lag_seconds <= 0:
        raise ValueError("maximum_display_lag_seconds must be positive.")
    if normal_fit_max_lag_seconds is not None and normal_fit_max_lag_seconds <= 0:
        raise ValueError("normal_fit_max_lag_seconds must be positive.")

    cell_msd = analysis["cell_msd"].copy()
    if msd_column not in cell_msd.columns:
        raise ValueError(f"MSD column {msd_column!r} is unavailable.")
    available_channels = sorted(cell_msd["channel"].dropna().astype(int).unique())
    selected_channels = _normalize_selection(available_channels, channels)
    labels = {
        int(key): str(value)
        for key, value in analysis["configuration"]["channel_condition_labels"].items()
    }

    calibration = (
        analysis["tracking_data"]
        .loc[
            analysis["tracking_data"]["channel"].isin(selected_channels),
            ["frame_interval_s", "is_3d"],
        ]
        .drop_duplicates()
    )
    frame_intervals = calibration["frame_interval_s"].dropna().unique()
    if len(frame_intervals) != 1:
        raise ValueError(
            "Paired channel MSD aggregation requires one common frame interval."
        )
    dimensions = calibration["is_3d"].dropna().astype(bool).unique()
    if len(dimensions) != 1:
        raise ValueError("Selected data mix 2D and 3D tracking.")
    dimension_factor = 6.0 if bool(dimensions[0]) else 4.0

    cell_ids = sorted(cell_msd["global_cell_id"].unique())
    if paired_cells_only:
        cell_ids = _complete_cell_ids(
            analysis["cell_summary"],
            selected_channels,
            require_valid_normal_fit=require_valid_normal_fit,
        )
    curves = cell_msd[
        cell_msd["global_cell_id"].isin(cell_ids)
        & cell_msd["channel"].isin(selected_channels)
        & cell_msd["minimum_contributors_met"].fillna(False)
    ].copy()
    curves[msd_column] = pd.to_numeric(curves[msd_column], errors="coerce")
    curves = curves[np.isfinite(curves[msd_column])]
    curve_keys = ["global_cell_id", "lag_frames", "lag_seconds", "channel"]
    if curves.duplicated(curve_keys).any():
        raise ValueError("Cell MSD table contains duplicate cell/channel/lag rows.")

    if paired_cells_only:
        wide = curves.pivot(
            index=["global_cell_id", "lag_frames", "lag_seconds"],
            columns="channel",
            values=msd_column,
        )
        for channel in selected_channels:
            if channel not in wide.columns:
                raise ValueError(f"Channel {channel} has no usable MSD curves.")
        wide = wide.dropna(subset=selected_channels, how="any")
        curves = (
            wide[selected_channels]
            .rename_axis(columns="channel")
            .stack()
            .rename(msd_column)
            .reset_index()
        )

    curves["channel_condition"] = curves["channel"].map(labels)
    rows: list[dict[str, Any]] = []
    for (channel, lag_frames, lag_seconds), group in curves.groupby(
        ["channel", "lag_frames", "lag_seconds"], sort=True
    ):
        values = group[msd_column].to_numpy(dtype=float)
        spread = _curve_spread(
            values,
            center=center,
            uncertainty=uncertainty,
        )
        rows.append(
            {
                "channel": int(channel),
                "channel_condition": labels[int(channel)],
                "lag_frames": int(lag_frames),
                "lag_seconds": float(lag_seconds),
                **spread,
            }
        )
    summary = pd.DataFrame(rows).sort_values(["channel", "lag_seconds"])
    if summary.empty:
        raise ValueError("No paired cell MSD values are available.")

    if normal_fit_max_lag_seconds is None:
        normal_fit_max_lag_seconds = analysis["configuration"].get("normal_fit_max_lag_seconds")
    if normal_fit_max_lag_seconds is None:
        normal_fit_max_lag_seconds = 10.0

    fit_rows: list[dict[str, Any]] = []
    if normal_fit_max_lag_seconds is not None:
        for channel in selected_channels:
            fit_data = summary[
                (summary["channel"] == channel)
                & (summary["lag_seconds"] > 0)
                & (summary["lag_seconds"] <= normal_fit_max_lag_seconds)
            ].dropna(subset=["center_msd_um2"])
            if len(fit_data) < 2:
                raise ValueError(
                    f"Channel {channel} has fewer than two group MSD points "
                    "in the requested fit interval."
                )
            regression = linregress(
                fit_data["lag_seconds"].to_numpy(dtype=float),
                fit_data["center_msd_um2"].to_numpy(dtype=float),
            )
            fit_rows.append(
                {
                    "channel": int(channel),
                    "channel_condition": labels[int(channel)],
                    "fit_points": int(len(fit_data)),
                    "fit_lag_start_s": float(fit_data["lag_seconds"].min()),
                    "fit_lag_end_s": float(fit_data["lag_seconds"].max()),
                    "dimension_factor": dimension_factor,
                    "D_group_curve_um2_s": float(regression.slope / dimension_factor),
                    "D_group_curve_se_um2_s": (
                        float(regression.stderr / dimension_factor)
                        if len(fit_data) > 2 and np.isfinite(regression.stderr)
                        else np.nan
                    ),
                    "slope_um2_s": float(regression.slope),
                    "offset_um2": float(regression.intercept),
                    "r_squared": float(regression.rvalue**2),
                }
            )
    fit_summary = pd.DataFrame(fit_rows)
    consistency_rows: list[dict[str, Any]] = []
    cell_d = analysis["cell_summary"].copy()
    cell_d = cell_d[
        cell_d["channel"].isin(selected_channels) & cell_d["D_fit_valid"].fillna(False)
    ]
    for channel in selected_channels:
        channel_d = cell_d[
            (cell_d["channel"] == channel) & cell_d["global_cell_id"].isin(cell_ids)
        ].copy()
        channel_curves = curves[curves["channel"] == channel]
        msd_cell_ids = set(channel_curves["global_cell_id"].unique())
        d_cell_ids = set(channel_d["global_cell_id"].unique())
        included_ids = sorted(d_cell_ids.intersection(msd_cell_ids))
        included_d = channel_d[channel_d["global_cell_id"].isin(included_ids)]
        group_fit = (
            fit_summary[fit_summary["channel"] == channel]
            if "channel" in fit_summary.columns
            else pd.DataFrame()
        )
        group_fit_d = (
            float(group_fit.iloc[0]["D_group_curve_um2_s"])
            if not group_fit.empty
            else np.nan
        )
        mean_cell_d = (
            float(included_d["D_um2_s"].mean()) if not included_d.empty else np.nan
        )
        median_cell_d = (
            float(included_d["D_um2_s"].median()) if not included_d.empty else np.nan
        )
        for cell_id in sorted(d_cell_ids.union(msd_cell_ids)):
            d_value = channel_d.loc[
                channel_d["global_cell_id"] == cell_id,
                "D_um2_s",
            ]
            consistency_rows.append(
                {
                    "global_cell_id": cell_id,
                    "channel": int(channel),
                    "channel_condition": labels[int(channel)],
                    "included_in_cell_D": cell_id in d_cell_ids,
                    "included_in_cell_MSD": cell_id in msd_cell_ids,
                    "included_in_both": cell_id in included_ids,
                    "cell_D_um2_s": (
                        float(d_value.iloc[0]) if not d_value.empty else np.nan
                    ),
                    "n_consistent_cells": int(len(included_ids)),
                    "mean_cell_D_um2_s": mean_cell_d,
                    "median_cell_D_um2_s": median_cell_d,
                    "D_group_curve_um2_s": group_fit_d,
                    "group_D_minus_mean_cell_D_um2_s": (
                        group_fit_d - mean_cell_d
                        if np.isfinite(group_fit_d) and np.isfinite(mean_cell_d)
                        else np.nan
                    ),
                    "identical_D_and_MSD_cell_sets": d_cell_ids == msd_cell_ids,
                }
            )
    consistency = pd.DataFrame(consistency_rows)
    plot_summary = summary.copy()
    plot_curves = curves.copy()
    if maximum_display_lag_seconds is not None:
        plot_summary = plot_summary[
            plot_summary["lag_seconds"] <= maximum_display_lag_seconds
        ]
        plot_curves = plot_curves[
            plot_curves["lag_seconds"] <= maximum_display_lag_seconds
        ]
    return {
        "cell_curves": curves,
        "summary": summary,
        "fit_summary": fit_summary,
        "D_msd_consistency": consistency,
        "plot_cell_curves": plot_curves,
        "plot_summary": plot_summary,
        "configuration": {
            "channels": selected_channels,
            "channel_condition_labels": labels,
            "msd_column": msd_column,
            "center": center,
            "uncertainty": uncertainty,
            "paired_cells_only": bool(paired_cells_only),
            "require_valid_normal_fit": bool(require_valid_normal_fit),
            "maximum_display_lag_seconds": maximum_display_lag_seconds,
            "normal_fit_max_lag_seconds": normal_fit_max_lag_seconds,
            "frame_interval_s": float(frame_intervals[0]),
            "dimension_factor": dimension_factor,
        },
    }


def _adjust_bh(p_values: Sequence[float]) -> np.ndarray:
    """Apply the Benjamini-Hochberg correction to finite p-values."""
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(values))
    if finite_indices.size == 0:
        return adjusted
    finite = values[finite_indices]
    order = np.argsort(finite)
    ranked = finite[order]
    n = len(ranked)
    scaled = ranked * n / np.arange(1, n + 1)
    monotonic = np.minimum.accumulate(scaled[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(monotonic, 0, 1)
    adjusted[finite_indices] = restored
    return adjusted


def compute_paired_channel_statistics(
    analysis: Mapping[str, Any],
    *,
    metrics: Sequence[str],
    channel_pair: Sequence[int] = (0, 1),
    test: str = "wilcoxon",
    log10_metrics: Sequence[str] = ("D_um2_s", "K_alpha_um2_s_alpha"),
    adjust_fdr: bool = True,
) -> pd.DataFrame:
    """Compare channel conditions using cells aligned by global_cell_id."""
    if len(channel_pair) != 2:
        raise ValueError("channel_pair must contain exactly two channels.")
    channel_a, channel_b = [int(value) for value in channel_pair]
    test = str(test).lower()
    if test not in {"wilcoxon", "paired_t", "paired_t_log10"}:
        raise ValueError("test must be 'wilcoxon', 'paired_t', or 'paired_t_log10'.")
    labels = analysis["configuration"]["channel_condition_labels"]
    data = analysis["cell_analysis_summary"]
    rows: list[dict[str, Any]] = []
    for metric in metrics:
        valid = _metric_valid_rows(data, metric)
        selected = valid[valid["channel"].isin([channel_a, channel_b])]
        pivot = selected.pivot(
            index="global_cell_id",
            columns="channel",
            values=metric,
        )
        if channel_a not in pivot.columns or channel_b not in pivot.columns:
            paired = pd.DataFrame(columns=[channel_a, channel_b])
        else:
            paired = pivot[[channel_a, channel_b]].dropna()
        raw_a = paired[channel_a].to_numpy(dtype=float)
        raw_b = paired[channel_b].to_numpy(dtype=float)
        transformed = False
        used_a = raw_a.copy()
        used_b = raw_b.copy()
        if test == "paired_t_log10" or (
            test == "paired_t" and metric in set(log10_metrics)
        ):
            positive = (used_a > 0) & (used_b > 0)
            used_a = np.log10(used_a[positive])
            used_b = np.log10(used_b[positive])
            raw_a = raw_a[positive]
            raw_b = raw_b[positive]
            transformed = True
        n = int(len(used_a))
        statistic = np.nan
        p_value = np.nan
        test_name = test
        if n >= 2:
            if test == "wilcoxon":
                differences = used_b - used_a
                if np.allclose(differences, 0):
                    statistic, p_value = 0.0, 1.0
                else:
                    result = wilcoxon(
                        used_a,
                        used_b,
                        alternative="two-sided",
                        zero_method="wilcox",
                    )
                    statistic = float(result.statistic)
                    p_value = float(result.pvalue)
            else:
                result = ttest_rel(used_a, used_b, nan_policy="omit")
                statistic = float(result.statistic)
                p_value = float(result.pvalue)
                test_name = "paired_t_log10" if transformed else "paired_t"
        ratios = (
            raw_b / raw_a
            if n and np.all(raw_a > 0) and np.all(raw_b > 0)
            else np.full(n, np.nan)
        )
        rows.append(
            {
                "metric": metric,
                "unit": METRIC_SPECS.get(metric, {}).get("unit", ""),
                "channel_1": channel_a,
                "condition_1": labels[channel_a],
                "channel_2": channel_b,
                "condition_2": labels[channel_b],
                "n_paired_cells": n,
                "test": test_name,
                "log10_transformed": transformed,
                "statistic": statistic,
                "p_value": p_value,
                "mean_condition_1": float(np.mean(raw_a)) if n else np.nan,
                "mean_condition_2": float(np.mean(raw_b)) if n else np.nan,
                "median_paired_difference_2_minus_1": (
                    float(np.median(raw_b - raw_a)) if n else np.nan
                ),
                "median_paired_ratio_2_over_1": (
                    float(np.nanmedian(ratios))
                    if n and np.isfinite(ratios).any()
                    else np.nan
                ),
            }
        )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["p_value_fdr_bh"] = (
            _adjust_bh(result["p_value"])
            if adjust_fdr
            else result["p_value"].to_numpy(dtype=float)
        )
    return result


def _metric_summary(
    data: pd.DataFrame,
    metric: str,
    channels: Sequence[int],
    *,
    center: str,
    uncertainty: str,
) -> pd.DataFrame:
    """Summarize finite cell-level values by channel and uncertainty choice."""
    rows: list[dict[str, Any]] = []
    for channel in channels:
        group = data[data["channel"] == channel]
        values = group[metric].to_numpy(dtype=float)
        n = int(values.size)
        if n == 0:
            continue
        mean = float(np.mean(values))
        median = float(np.median(values))
        center_value = mean if center == "mean" else median
        sd = float(np.std(values, ddof=1)) if n > 1 else np.nan
        sem = float(sd / math.sqrt(n)) if n > 1 else np.nan
        q1 = float(np.quantile(values, 0.25))
        q3 = float(np.quantile(values, 0.75))
        if uncertainty == "sem":
            lower, upper = center_value - sem, center_value + sem
        elif uncertainty == "sd":
            lower, upper = center_value - sd, center_value + sd
        elif uncertainty == "iqr":
            lower, upper = q1, q3
        else:
            lower = upper = center_value
        rows.append(
            {
                "metric": metric,
                "channel": int(channel),
                "channel_condition": group["channel_condition"].iloc[0],
                "n_cells": n,
                "mean": mean,
                "median": median,
                "center": center_value,
                "sd": sd,
                "sem": sem,
                "q1": q1,
                "q3": q3,
                "lower": lower,
                "upper": upper,
            }
        )
    return pd.DataFrame(rows)


def _format_p_value(value: float, *, symbol: str = "p") -> str:
    """Format a p-value for a plot annotation."""
    if not np.isfinite(value):
        return f"{symbol} = NA"
    if value < 0.0001:
        return f"{symbol} < 0.0001"
    return f"{symbol} = {value:.3g}"


def _significance_stars(value: float) -> str:
    """Convert a p-value to a conventional significance label."""
    if not np.isfinite(value):
        return "NA"
    if value < 0.0001:
        return "****"
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def _draw_significance_bracket(
    ax: plt.Axes,
    *,
    x1: float,
    x2: float,
    y_values: np.ndarray,
    label: str,
    log_scale: bool,
) -> None:
    """Draw a paired-comparison bracket above a metric plot."""
    finite = y_values[np.isfinite(y_values)]
    if finite.size == 0:
        return
    if log_scale:
        positive = finite[finite > 0]
        if positive.size == 0:
            return
        y = float(np.max(positive) * 1.22)
        h = y * 0.08
        ax.set_ylim(top=max(ax.get_ylim()[1], y + h * 2.6))
    else:
        y_min = float(np.min(finite))
        y_max = float(np.max(finite))
        span = max(y_max - y_min, abs(y_max) * 0.1, 1e-12)
        y = y_max + 0.13 * span
        h = 0.045 * span
        ax.set_ylim(top=max(ax.get_ylim()[1], y + h * 3.0))
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=1.3)
    ax.text(
        (x1 + x2) / 2,
        y + h * 1.2,
        label,
        ha="center",
        va="bottom",
        fontsize=mpl.rcParams["font.size"],
        fontweight="bold",
        color="black",
    )


def _draw_cell_metric_on_ax(
    ax: plt.Axes,
    analysis: Mapping[str, Any],
    *,
    metric: str,
    channels: Sequence[int],
    center: str,
    uncertainty: str,
    show_trajectory_background: bool,
    connect_paired_cells: bool,
    show_box: bool,
    box_width: float = 0.28,
    cap_width: float | None = 0.18,
    show_center_uncertainty: bool,
    condition_spacing: float,
    log_scale: bool,
    display_upper_percentile: float | None,
    statistics: pd.DataFrame | None,
    statistics_label_style: str,
    show_legend: bool,
    show_sample_counts: bool,
) -> pd.DataFrame:
    """Draw cell-level metric distributions and paired comparisons."""
    center, uncertainty = _validate_summary_options(center, uncertainty)
    if metric not in METRIC_SPECS:
        raise ValueError(f"Unsupported metric {metric!r}.")
    data = _metric_valid_rows(analysis["cell_analysis_summary"], metric)
    data = data[data["channel"].isin(channels)].copy()
    if log_scale:
        data = data[data[metric] > 0]
    if data.empty:
        raise ValueError(f"No valid cell-level values are available for {metric}.")
    if display_upper_percentile is not None and not (
        0 < float(display_upper_percentile) <= 100
    ):
        raise ValueError("display_upper_percentile must be in (0, 100].")
    if statistics_label_style not in {"full", "stars"}:
        raise ValueError("statistics_label_style must be 'full' or 'stars'.")

    labels = analysis["configuration"]["channel_condition_labels"]
    encodings = _channel_encoding(channels)
    if float(condition_spacing) <= 0:
        raise ValueError("condition_spacing must be positive.")
    positions = {
        channel: index * float(condition_spacing)
        for index, channel in enumerate(channels)
    }
    display_thresholds: dict[int, float] = {}
    display_parts: list[pd.DataFrame] = []
    for channel in channels:
        channel_data = data[data["channel"] == channel]
        values = channel_data[metric].to_numpy(dtype=float)
        threshold = (
            float(np.percentile(values, float(display_upper_percentile)))
            if display_upper_percentile is not None
            else np.inf
        )
        display_thresholds[channel] = threshold
        display_parts.append(channel_data[channel_data[metric] <= threshold])
    display_data = pd.concat(display_parts, ignore_index=False)
    all_cells = sorted(data["global_cell_id"].unique())
    jitter_rng = np.random.default_rng(417)
    cell_jitter = {
        cell: float(offset)
        for cell, offset in zip(
            all_cells,
            jitter_rng.uniform(-0.14, 0.14, max(len(all_cells), 1)),
        )
    }

    if show_trajectory_background:
        trajectory_metric = METRIC_SPECS[metric]["trajectory_metric"]
        trajectory = analysis.get("trajectory_spot_summary")
        if trajectory_metric and isinstance(trajectory, pd.DataFrame):
            for channel in channels:
                raw = trajectory.loc[
                    trajectory["channel"] == channel, trajectory_metric
                ]
                raw = pd.to_numeric(raw, errors="coerce")
                raw = raw[np.isfinite(raw)]
                if log_scale:
                    raw = raw[raw > 0]
                if raw.empty:
                    continue
                offsets = np.linspace(-0.22, 0.22, len(raw))
                ax.scatter(
                    positions[channel] + offsets,
                    raw,
                    s=14,
                    color="#A8A8A8",
                    alpha=0.20,
                    edgecolors="none",
                    zorder=0,
                )

    if connect_paired_cells and len(channels) >= 2:
        pivot = display_data.pivot(
            index="global_cell_id",
            columns="channel",
            values=metric,
        )
        complete = pivot.dropna(subset=channels, how="any")
        for cell_id, row in complete.iterrows():
            jitter = cell_jitter.get(cell_id, 0.0)
            xs = [positions[channel] + jitter for channel in channels]
            ys = [float(row[channel]) for channel in channels]
            ax.plot(
                xs,
                ys,
                color="#777777",
                alpha=0.34,
                linewidth=1.0,
                zorder=1,
            )

    if show_box:
        box_values = [
            display_data.loc[display_data["channel"] == channel, metric].to_numpy(
                dtype=float
            )
            for channel in channels
        ]
        bp_kwargs: dict[str, Any] = {
            "positions": [positions[channel] for channel in channels],
            "widths": float(box_width),
            "patch_artist": True,
            "showfliers": False,
            "whis": (5, 95),
            "zorder": 1,
        }
        if cap_width is not None:
            bp_kwargs["capwidths"] = float(cap_width)
        boxplot = ax.boxplot(box_values, **bp_kwargs)
        for patch in boxplot["boxes"]:
            patch.set_facecolor("white")
            patch.set_edgecolor("black")
            patch.set_linewidth(2.2)
        for element in ("whiskers", "caps", "medians"):
            for artist in boxplot[element]:
                artist.set_color("#D62728" if element == "medians" else "black")
                artist.set_linewidth(2.8 if element == "medians" else 2.2)

    for channel in channels:
        channel_data = display_data[display_data["channel"] == channel].sort_values(
            "global_cell_id"
        )
        x = np.asarray(
            [
                positions[channel] + cell_jitter.get(cell_id, 0.0)
                for cell_id in channel_data["global_cell_id"]
            ],
            dtype=float,
        )
        ax.scatter(
            x,
            channel_data[metric],
            s=34,
            marker="o",
            color="black",
            edgecolor="black",
            linewidth=0.35,
            alpha=0.95,
            zorder=3,
        )

    summary = _metric_summary(
        data,
        metric,
        channels,
        center=center,
        uncertainty=uncertainty,
    )
    summary["display_upper_percentile"] = display_upper_percentile
    summary["display_upper_threshold"] = summary["channel"].map(display_thresholds)
    summary["n_cells_above_display_percentile"] = summary["channel"].map(
        {
            channel: int(
                (
                    data.loc[data["channel"] == channel, metric]
                    > display_thresholds[channel]
                ).sum()
            )
            for channel in channels
        }
    )
    if show_center_uncertainty:
        for row in summary.itertuples(index=False):
            color, marker = encodings[int(row.channel)]
            lower_error = max(float(row.center - row.lower), 0.0)
            upper_error = max(float(row.upper - row.center), 0.0)
            yerr = (
                None
                if uncertainty == "none"
                else np.asarray([[lower_error], [upper_error]], dtype=float)
            )
            ax.errorbar(
                [positions[int(row.channel)]],
                [row.center],
                yerr=yerr,
                fmt=marker,
                markersize=12,
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=1.4,
                color="black",
                elinewidth=2.4,
                capsize=6,
                capthick=2.2,
                zorder=5,
            )

    if log_scale:
        ax.set_yscale("log")
    elif metric == "n_trajectories":
        ax.yaxis.set_major_locator(MaxNLocator(nbins=7, integer=True))
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=7, min_n_ticks=5))
    ax.set_xticks(
        [positions[channel] for channel in channels],
        [labels[channel] for channel in channels],
    )
    ax.set_ylabel(METRIC_SPECS[metric]["label"])
    ax.set_title(METRIC_SPECS[metric]["title"], pad=13)
    ax.set_xlabel("")
    ax.set_facecolor("white")
    ax.grid(False, which="both")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.6)
    ax.tick_params(axis="both", which="both", colors="black")
    font_family = mpl.rcParams["font.sans-serif"][0]
    ax.yaxis.label.set_fontfamily(font_family)
    ax.yaxis.label.set_color("black")
    for tick_label in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        tick_label.set_fontfamily(font_family)
        tick_label.set_color("black")
        tick_label.set_fontweight("semibold")
    if show_sample_counts:
        for row in summary.itertuples(index=False):
            hidden = int(row.n_cells_above_display_percentile)
            hidden_text = (
                f"; {hidden} >P{float(display_upper_percentile):g} hidden"
                if hidden and display_upper_percentile is not None
                else ""
            )
            ax.text(
                positions[int(row.channel)],
                -0.11,
                f"n = {int(row.n_cells)} cells{hidden_text}",
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=10.5,
                color="#555555",
            )

    if display_upper_percentile is not None:
        visible_upper = max(display_thresholds.values())
        if log_scale:
            ax.set_ylim(top=visible_upper * 1.18)
        else:
            visible_values = display_data[metric].to_numpy(dtype=float)
            span = max(
                visible_upper - float(np.min(visible_values)),
                abs(visible_upper) * 0.05,
                1e-12,
            )
            ax.set_ylim(top=visible_upper + span * 0.10)
    if (
        not log_scale
        and METRIC_SPECS[metric].get("zero_baseline", False)
        and float(display_data[metric].min()) >= 0
    ):
        ax.set_ylim(bottom=0)

    if statistics is not None and not statistics.empty and len(channels) >= 2:
        stat_row = statistics.loc[statistics["metric"] == metric]
        if not stat_row.empty:
            adjusted = "p_value_fdr_bh" in stat_row.columns
            p_column = "p_value_fdr_bh" if adjusted else "p_value"
            p_value = float(stat_row.iloc[0][p_column])
            significance = _significance_stars(p_value)
            label = (
                significance
                if statistics_label_style == "stars"
                else (
                    f"{significance}\n"
                    f"{_format_p_value(p_value, symbol='q' if adjusted else 'p')}; "
                    f"paired n = {int(stat_row.iloc[0]['n_paired_cells'])}"
                )
            )
            _draw_significance_bracket(
                ax,
                x1=positions[channels[0]],
                x2=positions[channels[1]],
                y_values=display_data[metric].to_numpy(dtype=float),
                label=label,
                log_scale=log_scale,
            )

    if show_legend:
        handles: list[Line2D] = [
            Line2D(
                [0],
                [0],
                marker=encodings[channel][1],
                color="none",
                markerfacecolor=encodings[channel][0],
                markeredgecolor="black",
                markersize=8,
                label=labels[channel],
            )
            for channel in channels
        ]
        if connect_paired_cells:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    color="#777777",
                    linewidth=1.2,
                    alpha=0.6,
                    label="Same cell",
                )
            )
        if show_trajectory_background and METRIC_SPECS[metric]["trajectory_metric"]:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="none",
                    markerfacecolor="#A8A8A8",
                    markeredgecolor="none",
                    alpha=0.45,
                    label="Trajectory (descriptive)",
                )
            )
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0,
            frameon=False,
        )
    return summary


def plot_cell_metric(
    analysis: Mapping[str, Any],
    *,
    metric: str,
    channels: Sequence[int] | int | None = None,
    center: str = "mean",
    uncertainty: str = "sem",
    show_trajectory_background: bool = False,
    connect_paired_cells: bool = False,
    show_box: bool = True,
    box_width: float = 0.32,
    cap_width: float | None = 0.20,
    show_center_uncertainty: bool = False,
    condition_spacing: float = 1.40,
    log_scale: bool | None = None,
    display_upper_percentile: float | None = 95.0,
    statistics: pd.DataFrame | None = None,
    statistics_label_style: str = "stars",
    show_legend: bool = False,
    show_sample_counts: bool = True,
    show_title: bool = True,
    ylabel: str | None = None,
    figsize: tuple[float, float] = (6.4, 6.2),
    plot_style: Mapping[str, Any] | None = None,
    output_directory: str | Path | None = None,
    show: bool = True,
    dpi: int = 600,
) -> tuple[plt.Figure, plt.Axes, pd.DataFrame]:
    """Plot a metric with cells as samples and channels as paired conditions."""
    center, uncertainty = _validate_summary_options(center, uncertainty)
    available = sorted(
        analysis["cell_analysis_summary"]["channel"].dropna().astype(int).unique()
    )
    selected_channels = _normalize_selection(available, channels)
    if log_scale is None:
        log_scale = bool(METRIC_SPECS[metric]["default_log"])
    with _publication_style(plot_style):
        fig, ax = plt.subplots(
            figsize=figsize,
            constrained_layout=True,
            facecolor="white",
        )
        summary = _draw_cell_metric_on_ax(
            ax,
            analysis,
            metric=metric,
            channels=selected_channels,
            center=center,
            uncertainty=uncertainty,
            show_trajectory_background=show_trajectory_background,
            connect_paired_cells=connect_paired_cells,
            show_box=show_box,
            box_width=box_width,
            cap_width=cap_width,
            show_center_uncertainty=show_center_uncertainty,
            condition_spacing=condition_spacing,
            log_scale=bool(log_scale),
            display_upper_percentile=display_upper_percentile,
            statistics=statistics,
            statistics_label_style=statistics_label_style,
            show_legend=show_legend,
            show_sample_counts=show_sample_counts,
        )
        if not show_title:
            ax.set_title("")
        if ylabel is not None:
            ax.set_ylabel(str(ylabel))
        if output_directory is not None:
            _save_figure(
                fig,
                output_directory=output_directory,
                stem=(
                    f"paired_{metric}_{center}_{uncertainty}_"
                    f"channels_{'_'.join(map(str, selected_channels))}"
                ),
                dpi=dpi,
            )
        if show:
            plt.show()
        else:
            plt.close(fig)
    return fig, ax, summary


def _draw_channel_msd_on_ax(
    ax: plt.Axes,
    msd_result: Mapping[str, Any],
    *,
    show_individual_cells: bool,
    error_style: str,
    scale: str,
    show_sample_counts: bool,
    single_label_legend: bool = True,
) -> tuple[list[Any], list[str]]:
    """Draw channel-level MSD curves, uncertainty, and optional fits."""
    summary = msd_result["plot_summary"]
    curves = msd_result["plot_cell_curves"]
    fit_summary = msd_result["fit_summary"]
    configuration = msd_result["configuration"]
    channels = [int(value) for value in configuration["channels"]]
    labels = configuration["channel_condition_labels"]
    encodings = _channel_encoding(channels)
    uncertainty = configuration["uncertainty"]
    center = configuration["center"]
    handles: list[Any] = []
    legend_labels: list[str] = []

    if show_individual_cells:
        for channel in channels:
            color, _ = encodings[channel]
            channel_curves = curves[curves["channel"] == channel]
            for _, cell_curve in channel_curves.groupby("global_cell_id"):
                ax.plot(
                    cell_curve["lag_seconds"],
                    cell_curve[configuration["msd_column"]],
                    color=color,
                    alpha=0.12,
                    linewidth=1.0,
                    zorder=0,
                )

    for channel in channels:
        color, marker = encodings[channel]
        data = summary[summary["channel"] == channel]
        if data.empty:
            continue
        lower = data["lower_msd_um2"].to_numpy(dtype=float)
        upper = data["upper_msd_um2"].to_numpy(dtype=float)
        center_values = data["center_msd_um2"].to_numpy(dtype=float)
        times = data["lag_seconds"].to_numpy(dtype=float)
        if error_style == "band" and uncertainty != "none":
            ax.fill_between(
                times,
                lower,
                upper,
                color=color,
                alpha=0.18,
                edgecolor="none",
                zorder=1,
            )
            (line,) = ax.plot(
                times,
                center_values,
                marker=marker,
                color=color,
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=1.1,
                markersize=8,
                linewidth=3.0,
                zorder=3,
            )
            mean_handle: Any = line
        else:
            yerr = None
            if uncertainty != "none":
                yerr = np.vstack(
                    [
                        np.maximum(center_values - lower, 0),
                        np.maximum(upper - center_values, 0),
                    ]
                )
            mean_handle = ax.errorbar(
                times,
                center_values,
                yerr=yerr,
                fmt=marker + "-",
                color=color,
                markerfacecolor=color,
                markeredgecolor="white",
                markeredgewidth=1.1,
                markersize=8,
                linewidth=3.0,
                elinewidth=1.8,
                capsize=4,
                capthick=1.7,
                zorder=3,
            )
        handles.append(mean_handle)

        fit = (
            fit_summary[fit_summary["channel"] == channel]
            if "channel" in fit_summary.columns
            else pd.DataFrame()
        )
        if not fit.empty:
            fit_row = fit.iloc[0]
            x_fit = np.linspace(
                float(fit_row["fit_lag_start_s"]),
                float(fit_row["fit_lag_end_s"]),
                200,
            )
            y_fit = float(fit_row["slope_um2_s"]) * x_fit + float(fit_row["offset_um2"])
            ax.plot(
                x_fit,
                y_fit,
                linestyle="--",
                color=color,
                linewidth=3.2,
                zorder=4,
            )
            d_val = float(fit_row["D_group_curve_um2_s"])
            if single_label_legend:
                legend_labels.append(f"{labels[channel]} (D={d_val:.2e} µm²/s)")
            else:
                summary_label = _summary_legend_label(center, uncertainty)
                if show_sample_counts:
                    n_min = int(data["n_cells"].min())
                    n_max = int(data["n_cells"].max())
                    n_text = f"n={n_min}" if n_min == n_max else f"n={n_min}–{n_max}"
                    legend_labels.append(
                        f"{labels[channel]}: {summary_label} ({n_text} cells)"
                    )
                else:
                    legend_labels.append(f"{labels[channel]}: {summary_label}")
                (fit_handle,) = ax.plot([], [], linestyle="--", color=color, linewidth=3.2)
                handles.append(fit_handle)
                legend_labels.append(
                    f"{labels[channel]} normal fit: D={d_val:.2e} µm²/s"
                )
        else:
            legend_labels.append(f"{labels[channel]}")

    ax.set_xlabel("Time lag (s)")
    ax.set_ylabel("MSD (µm²)")
    ax.set_title(f"{center.title()} MSD across paired cell samples", pad=13)
    ax.set_facecolor("white")
    ax.grid(False, which="both")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.6)
    ax.tick_params(axis="both", which="both", colors="black")
    if scale == "log":
        ax.set_xscale("log")
        ax.set_yscale("log")
    elif scale != "linear":
        raise ValueError("scale must be 'linear' or 'log'.")
    maximum_display = configuration["maximum_display_lag_seconds"]
    if maximum_display is not None and scale == "linear":
        ax.set_xlim(0, float(maximum_display) * 1.02)
    return handles, legend_labels


def plot_channel_msd(
    analysis: Mapping[str, Any],
    *,
    channels: Sequence[int] | int | None = None,
    msd_column: str = "pair_weighted_ensemble_msd_um2",
    center: str = "mean",
    uncertainty: str = "sem",
    paired_cells_only: bool = True,
    require_valid_normal_fit: bool = True,
    show_individual_cells: bool = False,
    error_style: str = "bars",
    scale: str = "linear",
    maximum_display_lag_seconds: float | None = None,
    normal_fit_max_lag_seconds: float | None = None,
    legend_outside_left: bool | None = None,
    legend_position: str | None = "top",
    figsize: tuple[float, float] | None = None,
    show_sample_counts: bool = True,
    show_title: bool = True,
    plot_style: Mapping[str, Any] | None = None,
    output_directory: str | Path | None = None,
    show: bool = True,
    dpi: int = 600,
) -> tuple[plt.Figure, plt.Axes, dict[str, Any]]:
    """Plot paired channel-condition MSD summaries across cell samples."""
    if error_style not in {"bars", "band"}:
        raise ValueError("error_style must be 'bars' or 'band'.")
    if legend_outside_left is not None:
        legend_position = "left" if legend_outside_left else "top"
    msd_result = build_channel_msd_summary(
        analysis,
        channels=channels,
        msd_column=msd_column,
        center=center,
        uncertainty=uncertainty,
        paired_cells_only=paired_cells_only,
        require_valid_normal_fit=require_valid_normal_fit,
        maximum_display_lag_seconds=maximum_display_lag_seconds,
        normal_fit_max_lag_seconds=normal_fit_max_lag_seconds,
    )
    with _publication_style(plot_style):
        if legend_position == "left":
            fig_size = figsize if figsize is not None else (13.5, 7.2)
            fig = plt.figure(
                figsize=fig_size,
                constrained_layout=True,
                facecolor="white",
            )
            layout = fig.add_gridspec(1, 2, width_ratios=(1.15, 2.35), wspace=0.06)
            legend_ax = fig.add_subplot(layout[0, 0])
            ax = fig.add_subplot(layout[0, 1])
            legend_ax.axis("off")
        elif legend_position == "top":
            fig_size = figsize if figsize is not None else (8.5, 6.8)
            fig, ax = plt.subplots(
                figsize=fig_size,
                constrained_layout=True,
                facecolor="white",
            )
            legend_ax = None
        else:
            fig_size = figsize if figsize is not None else (8.5, 6.5)
            fig, ax = plt.subplots(
                figsize=fig_size,
                constrained_layout=True,
                facecolor="white",
            )
            legend_ax = None
        handles, labels = _draw_channel_msd_on_ax(
            ax,
            msd_result,
            show_individual_cells=show_individual_cells,
            error_style=error_style,
            scale=scale,
            show_sample_counts=show_sample_counts,
        )
        if not show_title:
            ax.set_title("")
        leg_fs = (
            plot_style.get("legend.fontsize", mpl.rcParams.get("legend.fontsize", 10))
            if plot_style is not None
            else mpl.rcParams.get("legend.fontsize", 10)
        )
        if legend_position == "left" and legend_ax is not None:
            legend_ax.legend(
                handles,
                labels,
                loc="center left",
                frameon=False,
                fontsize=leg_fs,
                handlelength=2.7,
                labelspacing=0.9,
            )
        elif legend_position == "top":
            ax.legend(
                handles,
                labels,
                loc="lower center",
                bbox_to_anchor=(0.5, 1.02),
                ncol=len(handles) if len(handles) <= 3 else 2,
                frameon=False,
                fontsize=leg_fs,
                handlelength=2.5,
                columnspacing=1.5,
            )
        elif legend_position in {"right", "outside_right"}:
            ax.legend(
                handles,
                labels,
                loc="upper left",
                bbox_to_anchor=(1.02, 1.0),
                borderaxespad=0,
                frameon=False,
                fontsize=leg_fs,
            )
        elif legend_position is not None:
            ax.legend(
                handles,
                labels,
                loc=legend_position,
                frameon=False,
                fontsize=leg_fs,
            )
        if output_directory is not None:
            config = msd_result["configuration"]
            _save_figure(
                fig,
                output_directory=output_directory,
                stem=(
                    f"paired_channel_msd_{config['center']}_"
                    f"{config['uncertainty']}_{scale}"
                ),
                dpi=dpi,
            )
        if show:
            plt.show()
        else:
            plt.close(fig)
    return fig, ax, msd_result


def plot_combined_summary(
    analysis: Mapping[str, Any],
    *,
    panels: Sequence[str] = (
        "n_trajectories",
        "cell_spot_intensity_au",
        "cell_spot_fwhm_nm",
        "D_um2_s",
        "msd",
    ),
    channels: Sequence[int] | int | None = None,
    center: str = "mean",
    uncertainty: str = "sem",
    show_trajectory_background: bool = False,
    connect_paired_cells: bool = False,
    show_box: bool = True,
    box_width: float = 0.32,
    cap_width: float | None = 0.20,
    show_center_uncertainty: bool = False,
    condition_spacing: float = 1.40,
    display_upper_percentile: float | None = 95.0,
    msd_error_style: str = "band",
    msd_scale: str = "linear",
    maximum_display_lag_seconds: float | None = None,
    normal_fit_max_lag_seconds: float | None = None,
    single_label_legend: bool = True,
    statistics: pd.DataFrame | None = None,
    statistics_label_style: str = "stars",
    show_sample_counts: bool = True,
    show_panel_titles: bool = True,
    panel_ylabels: Mapping[str, str] | None = None,
    plot_style: Mapping[str, Any] | None = None,
    output_directory: str | Path | None = None,
    show: bool = True,
    dpi: int = 600,
) -> tuple[plt.Figure, list[plt.Axes]]:
    """Compose telomere-quality cell-metric panels and a paired MSD panel.

    Paired statistical brackets are drawn on cell-metric panels when
    ``statistics`` is provided. If omitted, the function uses
    ``analysis["paired_channel_statistics"]`` when available. Combined panels
    default to compact significance stars; exact adjusted q-values and paired
    sample sizes remain available in the statistics table.
    """
    available = sorted(
        analysis["cell_analysis_summary"]["channel"].dropna().astype(int).unique()
    )
    selected_channels = _normalize_selection(available, channels)
    invalid = [
        panel for panel in panels if panel != "msd" and panel not in METRIC_SPECS
    ]
    if invalid:
        raise ValueError(f"Unsupported combined-summary panels: {invalid}.")
    if not panels:
        raise ValueError("At least one panel must be requested.")
    if msd_error_style not in {"bars", "band"}:
        raise ValueError("msd_error_style must be 'bars' or 'band'.")
    if statistics is None:
        stored_statistics = analysis.get("paired_channel_statistics")
        if isinstance(stored_statistics, pd.DataFrame):
            statistics = stored_statistics
    msd_result = None
    if "msd" in panels:
        msd_result = build_channel_msd_summary(
            analysis,
            channels=selected_channels,
            center=center,
            uncertainty=uncertainty,
            maximum_display_lag_seconds=maximum_display_lag_seconds,
            normal_fit_max_lag_seconds=normal_fit_max_lag_seconds,
        )

    with _publication_style(plot_style):
        fig, axes_array = plt.subplots(
            1,
            len(panels),
            figsize=(5.8 * len(panels), 6.5),
            constrained_layout=True,
            facecolor="white",
        )
        axes = np.atleast_1d(axes_array).tolist()
        for ax, panel in zip(axes, panels):
            if panel == "msd":
                assert msd_result is not None
                handles, labels = _draw_channel_msd_on_ax(
                    ax,
                    msd_result,
                    show_individual_cells=False,
                    error_style=msd_error_style,
                    scale=msd_scale,
                    show_sample_counts=show_sample_counts,
                    single_label_legend=single_label_legend,
                )
                ax.legend(
                    handles,
                    labels,
                    loc="upper left",
                    fontsize=12,
                    frameon=True,
                    edgecolor="black",
                )
            else:
                _draw_cell_metric_on_ax(
                    ax,
                    analysis,
                    metric=panel,
                    channels=selected_channels,
                    center=center,
                    uncertainty=uncertainty,
                    show_trajectory_background=(
                        show_trajectory_background
                        and METRIC_SPECS[panel]["trajectory_metric"] is not None
                    ),
                    connect_paired_cells=connect_paired_cells,
                    show_box=show_box,
                    box_width=box_width,
                    cap_width=cap_width,
                    show_center_uncertainty=show_center_uncertainty,
                    condition_spacing=condition_spacing,
                    log_scale=bool(METRIC_SPECS[panel]["default_log"]),
                    display_upper_percentile=display_upper_percentile,
                    statistics=statistics,
                    statistics_label_style=statistics_label_style,
                    show_legend=False,
                    show_sample_counts=show_sample_counts,
                )
            if not show_panel_titles:
                ax.set_title("")
            if panel_ylabels is not None and panel in panel_ylabels:
                ax.set_ylabel(str(panel_ylabels[panel]))
        if "msd" in panels:
            # The MSD panel carries an x-axis label below its tick numbers while
            # the cell-metric panels end at their condition tick labels. Drop the
            # condition labels onto the same text row so every panel ends level.
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            msd_axis = axes[list(panels).index("msd")]
            target = msd_axis.xaxis.label.get_window_extent(renderer).y1
            for ax, panel in zip(axes, panels):
                if panel == "msd":
                    continue
                tick_labels = [
                    text for text in ax.get_xticklabels() if text.get_text()
                ]
                if not tick_labels:
                    continue
                top = max(
                    text.get_window_extent(renderer).y1 for text in tick_labels
                )
                ax.tick_params(
                    axis="x",
                    pad=(
                        ax.xaxis.majorTicks[0].get_pad()
                        + (top - target) * 72.0 / fig.dpi
                    ),
                )
        if output_directory is not None:
            _save_figure(
                fig,
                output_directory=output_directory,
                stem="paired_channel_combined_summary",
                dpi=dpi,
            )
        if show:
            plt.show()
        else:
            plt.close(fig)
    return fig, axes


__all__ = [
    "analyze_final_results",
    "build_channel_condition_summary",
    "build_channel_msd_summary",
    "compute_paired_channel_statistics",
    "fit_lag_seconds_from_frames",
    "plot_cell_metric",
    "plot_channel_msd",
    "plot_combined_summary",
    "save_final_tables",
    "summarize_spot_metrics",
]
