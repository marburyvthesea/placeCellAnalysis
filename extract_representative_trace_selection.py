from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal, Optional

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import zscore


@dataclass(frozen=True)
class FrequencySelectionResult:
    region_matches: pd.DataFrame
    cell_matches: pd.DataFrame
    cell_summary: pd.DataFrame
    cell_window_event_hz: pd.DataFrame
    cell_event_starts: pd.DataFrame


@dataclass(frozen=True)
class StrictTraceEventDetectionResult:
    peak_indices: np.ndarray
    onset_indices: np.ndarray
    offset_indices: np.ndarray
    smoothed_trace: np.ndarray


def _find_first_matching_column(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    lookup = {str(col).lower(): str(col) for col in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


def standardize_cell_labels(
    labels_df: Optional[pd.DataFrame],
    n_cells: int,
    *,
    labels_are_one_based: bool = True,
    cell_id_col: Optional[str] = None,
) -> pd.DataFrame:
    if labels_df is None:
        return pd.DataFrame(
            {
                "cell_id": np.arange(1, n_cells + 1, dtype=int),
                "cell_index": np.arange(n_cells, dtype=int),
            }
        )

    standardized = labels_df.copy()

    if "cell_index" not in standardized.columns:
        if cell_id_col is None:
            cell_id_col = _find_first_matching_column(
                standardized.columns,
                ["cell_id", "cell", "roi", "roi_id", "component", "unit", "neuron", "index"],
            )
        if cell_id_col is None:
            standardized["cell_id"] = np.arange(1, len(standardized) + 1, dtype=int)
            cell_id_col = "cell_id"
        standardized["cell_id"] = standardized[cell_id_col].astype(int)
        if labels_are_one_based:
            standardized["cell_index"] = standardized["cell_id"] - 1
        else:
            standardized["cell_index"] = standardized["cell_id"]
    else:
        standardized["cell_index"] = standardized["cell_index"].astype(int)
        if "cell_id" not in standardized.columns:
            if labels_are_one_based:
                standardized["cell_id"] = standardized["cell_index"] + 1
            else:
                standardized["cell_id"] = standardized["cell_index"]

    standardized = standardized[(standardized["cell_index"] >= 0) & (standardized["cell_index"] < n_cells)].copy()
    standardized = standardized.drop_duplicates(subset=["cell_index"]).sort_values("cell_index").reset_index(drop=True)
    standardized["cell_id"] = standardized["cell_id"].astype(int)
    return standardized


def zscore_traces_by_cell(traces_2d: np.ndarray) -> np.ndarray:
    traces_2d = np.asarray(traces_2d, dtype=np.float32)
    z_traces = zscore(traces_2d, axis=0, nan_policy="omit")
    z_traces = np.asarray(z_traces, dtype=np.float32)
    z_traces[~np.isfinite(z_traces)] = 0.0
    return z_traces


def _smooth_trace_centered(trace_1d: np.ndarray, window_samples: Optional[int]) -> np.ndarray:
    trace_1d = np.asarray(trace_1d, dtype=np.float32).reshape(-1)
    if window_samples is None or int(window_samples) <= 1:
        return trace_1d.copy()
    return (
        pd.Series(trace_1d)
        .rolling(window=int(window_samples), center=True, min_periods=1)
        .mean()
        .to_numpy(dtype=np.float32)
    )


def _iter_with_optional_progress(
    iterable,
    *,
    show_progress: bool = False,
    total: Optional[int] = None,
    desc: Optional[str] = None,
    unit: str = "item",
    leave: bool = True,
):
    if not show_progress:
        return iterable

    try:
        from tqdm.auto import tqdm
    except ImportError:
        return iterable

    return tqdm(
        iterable,
        total=total,
        desc=desc,
        unit=unit,
        leave=leave,
        dynamic_ncols=True,
    )


def compute_component_centroids(spatial_3d: np.ndarray, cell_indices: np.ndarray) -> np.ndarray:
    spatial_3d = np.asarray(spatial_3d, dtype=np.float32)
    if spatial_3d.ndim != 3:
        raise ValueError(f"Expected spatial_3d to be 3D, got shape {spatial_3d.shape}")

    image_h, image_w, n_cells = spatial_3d.shape
    invalid = cell_indices[(cell_indices < 0) | (cell_indices >= n_cells)]
    if invalid.size:
        raise ValueError(f"Some cell indices are out of bounds for spatial_3d: {invalid[:10]}")

    ys, xs = np.mgrid[0:image_h, 0:image_w]
    centroids = np.full((len(cell_indices), 2), np.nan, dtype=np.float32)

    for out_idx, cell_index in enumerate(cell_indices):
        component = spatial_3d[:, :, cell_index].copy()
        component[component < 0] = 0
        total = component.sum()
        if total <= 0:
            continue
        centroids[out_idx, 0] = float((xs * component).sum() / total)
        centroids[out_idx, 1] = float((ys * component).sum() / total)

    return centroids


def assign_grid_regions(
    centroids_xy: np.ndarray,
    image_shape: tuple[int, int],
    region_grid_shape: tuple[int, int] = (4, 4),
) -> pd.DataFrame:
    image_h, image_w = image_shape
    n_rows, n_cols = region_grid_shape
    if n_rows < 1 or n_cols < 1:
        raise ValueError("region_grid_shape entries must both be >= 1")

    x_edges = np.linspace(0, image_w, n_cols + 1)
    y_edges = np.linspace(0, image_h, n_rows + 1)

    region_rows = np.full(len(centroids_xy), -1, dtype=int)
    region_cols = np.full(len(centroids_xy), -1, dtype=int)
    region_ids = np.full(len(centroids_xy), "unassigned", dtype=object)

    for idx, (x_coord, y_coord) in enumerate(centroids_xy):
        if not np.isfinite(x_coord) or not np.isfinite(y_coord):
            continue

        region_col = int(np.clip(np.searchsorted(x_edges, x_coord, side="right") - 1, 0, n_cols - 1))
        region_row = int(np.clip(np.searchsorted(y_edges, y_coord, side="right") - 1, 0, n_rows - 1))

        region_rows[idx] = region_row
        region_cols[idx] = region_col
        region_ids[idx] = f"r{region_row}_c{region_col}"

    return pd.DataFrame(
        {
            "centroid_x": centroids_xy[:, 0],
            "centroid_y": centroids_xy[:, 1],
            "region_row": region_rows,
            "region_col": region_cols,
            "region_id": region_ids,
        }
    )


def detect_event_starts_from_z_trace(
    z_trace: np.ndarray,
    *,
    event_threshold_z: float = 2.5,
    min_peak_distance_frames: int = 20,
    min_peak_prominence: Optional[float] = None,
) -> np.ndarray:
    z_trace = np.asarray(z_trace, dtype=np.float32).reshape(-1)
    peak_kwargs: dict[str, float | int] = {
        "height": float(event_threshold_z),
        "distance": max(1, int(min_peak_distance_frames)),
    }
    if min_peak_prominence is not None:
        peak_kwargs["prominence"] = float(min_peak_prominence)

    peak_indices, _ = find_peaks(z_trace, **peak_kwargs)
    return peak_indices.astype(int, copy=False)


def detect_strict_peak_stats_events_from_z_trace(
    z_trace: np.ndarray,
    *,
    threshold_z: float = 2.5,
    smooth_window_samples: int = 3,
    min_peak_prominence: float = 1.5,
    min_peak_distance_frames: int = 20,
    min_peak_width_frames: int = 2,
    return_to_baseline_tol: float = 0.4,
    return_to_baseline_window_frames: int = 5,
    onset_tol: float = 0.25,
    onset_window_len_frames: int = 3,
    offset_tol: float = 0.25,
    offset_window_len_frames: int = 3,
) -> StrictTraceEventDetectionResult:
    """Replicate the stricter peak-stats logic on one z-scored trace."""
    z_trace = np.asarray(z_trace, dtype=np.float32).reshape(-1)
    n_frames = z_trace.size
    if n_frames == 0:
        empty = np.array([], dtype=int)
        return StrictTraceEventDetectionResult(
            peak_indices=empty,
            onset_indices=empty,
            offset_indices=empty,
            smoothed_trace=z_trace.copy(),
        )

    smoothed_trace = _smooth_trace_centered(z_trace, smooth_window_samples)
    peak_indices, _ = find_peaks(
        smoothed_trace,
        height=float(threshold_z),
        distance=max(1, int(min_peak_distance_frames)),
        prominence=float(min_peak_prominence),
        width=float(min_peak_width_frames),
    )
    peak_indices = peak_indices.astype(int, copy=False)

    if peak_indices.size <= 1:
        accepted_peak_indices = peak_indices
    else:
        baseline_window = max(1, int(return_to_baseline_window_frames))
        accepted = [int(peak_indices[0])]

        def has_baseline_return(start_idx: int, stop_idx: int) -> bool:
            if stop_idx <= start_idx + 1:
                return True
            segment = np.abs(smoothed_trace[start_idx + 1:stop_idx])
            if segment.size == 0:
                return True
            if segment.size < baseline_window:
                return bool(np.all(segment <= float(return_to_baseline_tol)))

            baseline_mask = segment <= float(return_to_baseline_tol)
            kernel = np.ones(baseline_window, dtype=int)
            return bool(np.any(np.convolve(baseline_mask.astype(int), kernel, mode="valid") == baseline_window))

        for candidate_idx in peak_indices[1:]:
            if has_baseline_return(accepted[-1], int(candidate_idx)):
                accepted.append(int(candidate_idx))

        accepted_peak_indices = np.asarray(accepted, dtype=int)

    onset_indices = np.zeros(accepted_peak_indices.size, dtype=int)
    offset_indices = np.full(accepted_peak_indices.size, n_frames - 1, dtype=int)
    offset_window_len_frames = max(1, int(offset_window_len_frames))
    onset_window_len_frames = max(1, int(onset_window_len_frames))
    last_offset_start = max(0, n_frames - offset_window_len_frames)

    for event_idx, peak_idx in enumerate(accepted_peak_indices):
        onset_idx = 0
        for trace_idx in range(int(peak_idx), onset_window_len_frames - 2, -1):
            start_idx = trace_idx - onset_window_len_frames + 1
            if start_idx < 0:
                break
            if np.all(np.abs(z_trace[start_idx:trace_idx + 1]) < float(onset_tol)):
                onset_idx = start_idx
                break

        baseline_value = float(z_trace[onset_idx])
        offset_idx = n_frames - 1
        for trace_idx in range(int(peak_idx), last_offset_start + 1):
            window = slice(trace_idx, trace_idx + offset_window_len_frames)
            if np.all(np.abs(z_trace[window]) < float(offset_tol)) and np.all(z_trace[window] < baseline_value):
                offset_idx = trace_idx + offset_window_len_frames - 1
                break

        onset_indices[event_idx] = int(onset_idx)
        offset_indices[event_idx] = int(offset_idx)

    return StrictTraceEventDetectionResult(
        peak_indices=accepted_peak_indices,
        onset_indices=onset_indices,
        offset_indices=offset_indices,
        smoothed_trace=smoothed_trace,
    )


def build_event_start_matrix(
    traces_2d: np.ndarray,
    *,
    zscore_each_trace: bool = True,
    event_detection_method: Literal["simple", "strict_peak_stats"] = "simple",
    event_index_kind: Literal["peak", "onset"] = "peak",
    event_threshold_z: float = 2.5,
    min_peak_distance_frames: int = 20,
    min_peak_prominence: Optional[float] = None,
    strict_smooth_window_samples: int = 3,
    strict_min_peak_prominence: float = 1.5,
    strict_min_peak_width_frames: int = 2,
    strict_return_to_baseline_tol: float = 0.4,
    strict_return_to_baseline_window_frames: int = 5,
    strict_onset_tol: float = 0.25,
    strict_onset_window_len_frames: int = 3,
    strict_offset_tol: float = 0.25,
    strict_offset_window_len_frames: int = 3,
    show_progress: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    traces_2d = np.asarray(traces_2d, dtype=np.float32)
    if traces_2d.ndim != 2:
        raise ValueError(f"Expected traces_2d to be 2D, got shape {traces_2d.shape}")

    work_traces = zscore_traces_by_cell(traces_2d) if zscore_each_trace else traces_2d.copy()
    n_frames, n_cells = work_traces.shape
    event_start_matrix = np.zeros((n_frames, n_cells), dtype=np.uint8)

    cell_indices = _iter_with_optional_progress(
        range(n_cells),
        show_progress=show_progress,
        total=n_cells,
        desc="Detecting cell events",
        unit="cell",
        leave=True,
    )
    for cell_idx in cell_indices:
        if event_detection_method == "strict_peak_stats":
            strict_events = detect_strict_peak_stats_events_from_z_trace(
                work_traces[:, cell_idx],
                threshold_z=event_threshold_z,
                smooth_window_samples=strict_smooth_window_samples,
                min_peak_prominence=strict_min_peak_prominence,
                min_peak_distance_frames=min_peak_distance_frames,
                min_peak_width_frames=strict_min_peak_width_frames,
                return_to_baseline_tol=strict_return_to_baseline_tol,
                return_to_baseline_window_frames=strict_return_to_baseline_window_frames,
                onset_tol=strict_onset_tol,
                onset_window_len_frames=strict_onset_window_len_frames,
                offset_tol=strict_offset_tol,
                offset_window_len_frames=strict_offset_window_len_frames,
            )
            if event_index_kind == "onset":
                event_indices = strict_events.onset_indices
            else:
                event_indices = strict_events.peak_indices
        else:
            event_indices = detect_event_starts_from_z_trace(
                work_traces[:, cell_idx],
                event_threshold_z=event_threshold_z,
                min_peak_distance_frames=min_peak_distance_frames,
                min_peak_prominence=min_peak_prominence,
            )

        event_start_matrix[event_indices, cell_idx] = 1

    return work_traces, event_start_matrix


def compute_windowed_event_frequency(
    event_start_matrix: np.ndarray,
    *,
    fps: float = 20.0,
    window_size_frames: int = 1000,
    window_step_frames: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    event_start_matrix = np.asarray(event_start_matrix, dtype=np.uint8)
    if event_start_matrix.ndim != 2:
        raise ValueError(f"Expected event_start_matrix to be 2D, got shape {event_start_matrix.shape}")

    n_frames, _ = event_start_matrix.shape
    if window_size_frames < 1 or window_size_frames > n_frames:
        raise ValueError(
            f"window_size_frames must be between 1 and the number of frames ({n_frames}), got {window_size_frames}"
        )
    if window_step_frames < 1:
        raise ValueError("window_step_frames must be >= 1")

    cumulative = np.vstack(
        [
            np.zeros((1, event_start_matrix.shape[1]), dtype=np.int64),
            np.cumsum(event_start_matrix.astype(np.int64), axis=0),
        ]
    )
    full_window_counts = cumulative[window_size_frames:] - cumulative[:-window_size_frames]
    full_window_starts = np.arange(full_window_counts.shape[0], dtype=int)

    sampled_window_starts = full_window_starts[::window_step_frames]
    sampled_window_counts = full_window_counts[::window_step_frames]
    sampled_window_hz = sampled_window_counts / (window_size_frames / float(fps))
    return sampled_window_starts, sampled_window_hz.astype(np.float32)


def find_representative_cells_and_regions_by_frequency(
    traces_2d: np.ndarray,
    spatial_3d: np.ndarray,
    *,
    labels_df: Optional[pd.DataFrame] = None,
    target_event_hz: float,
    fps: float = 20.0,
    window_size_frames: int = 1000,
    window_step_frames: int = 50,
    zscore_each_trace: bool = True,
    event_detection_method: Literal["simple", "strict_peak_stats"] = "simple",
    event_index_kind: Literal["peak", "onset"] = "peak",
    event_threshold_z: float = 2.5,
    min_peak_distance_frames: int = 20,
    min_peak_prominence: Optional[float] = None,
    strict_smooth_window_samples: int = 3,
    strict_min_peak_prominence: float = 1.5,
    strict_min_peak_width_frames: int = 2,
    strict_return_to_baseline_tol: float = 0.4,
    strict_return_to_baseline_window_frames: int = 5,
    strict_onset_tol: float = 0.25,
    strict_onset_window_len_frames: int = 3,
    strict_offset_tol: float = 0.25,
    strict_offset_window_len_frames: int = 3,
    region_grid_shape: tuple[int, int] = (4, 4),
    min_cells_per_region: int = 3,
    top_k_regions: int = 5,
    top_k_cells_per_region: int = 4,
    labels_are_one_based: bool = True,
    cell_id_col: Optional[str] = None,
    show_progress: bool = False,
) -> FrequencySelectionResult:
    traces_2d = np.asarray(traces_2d, dtype=np.float32)
    if traces_2d.ndim != 2:
        raise ValueError(f"Expected traces_2d to be 2D, got shape {traces_2d.shape}")

    n_frames, n_cells_total = traces_2d.shape
    cell_table = standardize_cell_labels(
        labels_df,
        n_cells_total,
        labels_are_one_based=labels_are_one_based,
        cell_id_col=cell_id_col,
    ).copy()

    if cell_table.empty:
        raise ValueError("No valid cells remain after standardizing labels_df")

    selected_cell_indices = cell_table["cell_index"].to_numpy(dtype=int)
    selected_traces = traces_2d[:, selected_cell_indices]

    _, event_start_matrix = build_event_start_matrix(
        selected_traces,
        zscore_each_trace=zscore_each_trace,
        event_detection_method=event_detection_method,
        event_index_kind=event_index_kind,
        event_threshold_z=event_threshold_z,
        min_peak_distance_frames=min_peak_distance_frames,
        min_peak_prominence=min_peak_prominence,
        strict_smooth_window_samples=strict_smooth_window_samples,
        strict_min_peak_prominence=strict_min_peak_prominence,
        strict_min_peak_width_frames=strict_min_peak_width_frames,
        strict_return_to_baseline_tol=strict_return_to_baseline_tol,
        strict_return_to_baseline_window_frames=strict_return_to_baseline_window_frames,
        strict_onset_tol=strict_onset_tol,
        strict_onset_window_len_frames=strict_onset_window_len_frames,
        strict_offset_tol=strict_offset_tol,
        strict_offset_window_len_frames=strict_offset_window_len_frames,
        show_progress=show_progress,
    )

    window_starts, cell_window_event_hz = compute_windowed_event_frequency(
        event_start_matrix,
        fps=fps,
        window_size_frames=window_size_frames,
        window_step_frames=window_step_frames,
    )

    centroids_xy = compute_component_centroids(spatial_3d, selected_cell_indices)
    region_df = assign_grid_regions(centroids_xy, spatial_3d.shape[:2], region_grid_shape=region_grid_shape)

    cell_table = pd.concat([cell_table.reset_index(drop=True), region_df.reset_index(drop=True)], axis=1)

    abs_error_matrix = np.abs(cell_window_event_hz - float(target_event_hz))
    best_window_idx_per_cell = np.argmin(abs_error_matrix, axis=0)
    best_window_start_per_cell = window_starts[best_window_idx_per_cell]
    best_event_hz_per_cell = cell_window_event_hz[best_window_idx_per_cell, np.arange(cell_window_event_hz.shape[1])]
    best_abs_error_per_cell = abs_error_matrix[best_window_idx_per_cell, np.arange(abs_error_matrix.shape[1])]

    cell_summary = cell_table.copy()
    cell_summary["best_window_index"] = best_window_idx_per_cell
    cell_summary["best_window_start"] = best_window_start_per_cell
    cell_summary["best_window_stop_exclusive"] = best_window_start_per_cell + int(window_size_frames)
    cell_summary["best_event_hz"] = best_event_hz_per_cell
    cell_summary["best_abs_error_hz"] = best_abs_error_per_cell
    cell_summary = cell_summary.sort_values(["best_abs_error_hz", "best_event_hz"]).reset_index(drop=True)

    region_records: list[dict[str, object]] = []
    grouped_regions = list(cell_table.groupby("region_id", sort=False))
    grouped_regions = [
        (region_id, region_cells)
        for region_id, region_cells in grouped_regions
        if region_id != "unassigned" and len(region_cells) >= min_cells_per_region
    ]
    grouped_regions_iter = _iter_with_optional_progress(
        grouped_regions,
        show_progress=show_progress,
        total=len(grouped_regions),
        desc="Scoring spatial regions",
        unit="region",
        leave=False,
    )
    for region_id, region_cells in grouped_regions_iter:
        region_local_indices = region_cells.index.to_numpy(dtype=int)
        region_window_hz = cell_window_event_hz[:, region_local_indices]
        region_mean_hz = np.nanmean(region_window_hz, axis=1)
        region_mean_abs_error = np.abs(region_mean_hz - float(target_event_hz))

        best_region_window_idx = int(np.argmin(region_mean_abs_error))
        best_region_window_start = int(window_starts[best_region_window_idx])

        region_records.append(
            {
                "region_id": region_id,
                "region_row": int(region_cells["region_row"].iloc[0]),
                "region_col": int(region_cells["region_col"].iloc[0]),
                "window_index": best_region_window_idx,
                "window_start": best_region_window_start,
                "window_stop_exclusive": best_region_window_start + int(window_size_frames),
                "window_duration_seconds": float(window_size_frames / fps),
                "mean_event_hz": float(region_mean_hz[best_region_window_idx]),
                "mean_abs_error_hz": float(region_mean_abs_error[best_region_window_idx]),
                "n_region_cells": int(len(region_cells)),
                "region_cell_ids": region_cells["cell_id"].astype(int).tolist(),
                "region_cell_indices": region_cells["cell_index"].astype(int).tolist(),
            }
        )

    region_matches = pd.DataFrame(region_records)
    if region_matches.empty:
        raise ValueError(
            "No populated regions met min_cells_per_region. Try reducing min_cells_per_region or region_grid_shape."
        )

    region_matches = region_matches.sort_values(["mean_abs_error_hz", "n_region_cells"], ascending=[True, False]).reset_index(drop=True)
    region_matches.insert(0, "region_rank", np.arange(1, len(region_matches) + 1))
    region_matches = region_matches.head(top_k_regions).copy()

    cell_match_records: list[dict[str, object]] = []
    region_match_rows = list(region_matches.iterrows())
    region_match_rows_iter = _iter_with_optional_progress(
        region_match_rows,
        show_progress=show_progress,
        total=len(region_match_rows),
        desc="Selecting representative cells",
        unit="region",
        leave=False,
    )
    for _, region_row in region_match_rows_iter:
        region_id = str(region_row["region_id"])
        window_index = int(region_row["window_index"])
        region_cells = cell_table[cell_table["region_id"] == region_id].copy()
        region_local_indices = region_cells.index.to_numpy(dtype=int)
        region_event_hz = cell_window_event_hz[window_index, region_local_indices]
        region_abs_error = np.abs(region_event_hz - float(target_event_hz))

        region_cells["window_index"] = window_index
        region_cells["window_start"] = int(region_row["window_start"])
        region_cells["window_stop_exclusive"] = int(region_row["window_stop_exclusive"])
        region_cells["window_duration_seconds"] = float(region_row["window_duration_seconds"])
        region_cells["event_hz"] = region_event_hz
        region_cells["abs_error_hz"] = region_abs_error
        region_cells["region_rank"] = int(region_row["region_rank"])

        region_cells = region_cells.sort_values(["abs_error_hz", "event_hz"]).head(top_k_cells_per_region)
        region_cells = region_cells.reset_index(drop=True)
        region_cells.insert(len(region_cells.columns), "cell_rank_within_region", np.arange(1, len(region_cells) + 1))

        cell_match_records.extend(region_cells.to_dict(orient="records"))

    cell_matches = pd.DataFrame(cell_match_records)
    if not cell_matches.empty:
        cell_matches = cell_matches.sort_values(["region_rank", "cell_rank_within_region"]).reset_index(drop=True)

    cell_window_event_hz_df = pd.DataFrame(
        cell_window_event_hz,
        index=pd.Index(window_starts, name="window_start"),
        columns=cell_table["cell_id"].astype(int).tolist(),
    )
    event_start_df = pd.DataFrame(
        event_start_matrix,
        columns=cell_table["cell_id"].astype(int).tolist(),
    )

    return FrequencySelectionResult(
        region_matches=region_matches,
        cell_matches=cell_matches,
        cell_summary=cell_summary,
        cell_window_event_hz=cell_window_event_hz_df,
        cell_event_starts=event_start_df,
    )
