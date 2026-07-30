"""
FRAP Plotting Module for UTag Paper (Figure 4)
===============================================

Publication-quality plotting functions for FRAP analysis. Extracted from
FRAP_representative_images.ipynb and Figure_4_FRAP_interpreatation.ipynb
to eliminate code duplication and improve notebook readability.

Functions (Representative Images):
    - get_data_folder_path
    - plot_images_frap_all_channels_representative
    - plot_kymograph
    - plot_kymograph_downsampled
    - plot_merged_image
    - save_video_as_avi
    - compose_pngs

Functions (Interpretation / Statistics):
    - downsample_to_5_seconds
    - plot_FRAP_trajectories
    - save_t_half_source_data
    - write_frap_source_data_readme
    - plot_t_half_box_swarm
    - plot_mean_trajectories_all
    - plot_box_swarm_final_values

Usage:
    from frap_plotting_module import (
        plot_images_frap_all_channels_representative,
        plot_kymograph,
        plot_FRAP_trajectories,
        plot_t_half_box_swarm,
    )

Author: Luis Aguilera, Ning Zhao
"""

# Standard library
import io
from itertools import combinations
# Scientific computing
import numpy as np
import pandas as pd
import cv2
from scipy.optimize import curve_fit
import scipy.stats as stats
from scipy.ndimage import gaussian_filter

# Visualization
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib_scalebar.scalebar import ScaleBar
from skimage import exposure
from PIL import Image
import seaborn as sns


# ============================================================================
# Representative Images Functions
# ============================================================================

def get_data_folder_path(dataset, path_main_folder):
    """Return (data_folder_path, selected_image_name, list_axis_limits, y_label_list, plot_name) for a dataset."""
    configs = {
        'utag': {
            'subfolder': "FRAP 1xUtag",
            'lif': "20231207_FRAP_pNZ112_pNZ122.lif",
            'image': 'Cell05',
            'start_x': 150, 'end_y': 110,
            'ch0': 'Anti-UTag-FB-mEGFP', 'ch1': '1xUTag-mCh-H2B',
            'plot_name': 'Utag',
        },
        'utag_cf': {
            'subfolder': "FRAP 1xUTag_CF",
            'lif': "20250507 pNZ122_244 FRAP.lif",
            'image': 'FRAP 026',
            'start_x': 150, 'end_y': 100,
            'ch0': r'Anti-UTag-FB($\Delta$Cys)-mEGFP', 'ch1': '1xUTag-mCh-H2B',
            'plot_name': 'Utag_cf',
        },
        'suntag': {
            'subfolder': "FRAP 1xSUN",
            'lif': "20240306_pNZ317_219_FRAP.lif",
            'image': 'FRAP 015',
            'start_x': 150, 'end_y': 100,
            'ch0': 'Anti-SunTag-scFv-mEGFP', 'ch1': '1xSunTag-mCh-H2B',
            'plot_name': 'SunTag',
        },
        'alfatag': {
            'subfolder': "FRAP 1xALFA",
            'lif': "20240731_pNZ165_218.lif",
            'image': 'FRAP 011',
            'start_x': 130, 'end_y': 150,
            'ch0': 'Anti-ALFA-tag-nb-mEGFP', 'ch1': '1xALFA-tag-mCh-H2B',
            'plot_name': 'AlfaTag',
        },
        'hatag': {
            'subfolder': "FRAP 1xHA",
            'lif': "20240313_pNZ035_043_FRAP.lif",
            'image': 'FRAP 004',
            'start_x': 150, 'end_y': 100,
            'ch0': 'Anti-HA-FB-mEGFP', 'ch1': '1xHA-mCh-H2B',
            'plot_name': 'HA',
        },
    }
    cfg = configs[dataset]
    data_folder_path = path_main_folder / cfg['subfolder'] / cfg['lif']
    start_x = cfg['start_x']; end_x = start_x + 250
    end_y = cfg['end_y']; start_y = end_y + 250
    list_axis_limits = [start_x, end_x, start_y, end_y]
    y_label_list = [cfg['ch0'], cfg['ch1']]
    return data_folder_path, cfg['image'], list_axis_limits, y_label_list, cfg['plot_name']


def plot_images_frap_all_channels_representative(
    image_TZXYC, results_folder, pixel_xy_um, scalebar_size=0,
    list_selected_frames=[0, 10, 40, 100, 139],
    cmap_list=None, selected_color_channel=None, coordinates_roi=None,
    radius_roi_size_px=10, plot_name='temp.png', list_axis_limits=None,
    y_label_list=None, list_selected_frame_values_real_time=None,
    x_title_list=None, masks_TXY=None, min_size_image=150
):
    """Plot representative frames for all color channels with optional mask cropping."""
    if selected_color_channel is not None:
        image_TZXYC = image_TZXYC[..., selected_color_channel]
        image_TZXYC = np.expand_dims(image_TZXYC, axis=-1)

    number_color_channels = image_TZXYC.shape[-1]
    num_frames = len(list_selected_frames)

    # Crop window from mask
    if masks_TXY is not None:
        mask_xy = masks_TXY.max(axis=0)
        ys, xs = np.where(mask_xy > 0)
        if ys.size == 0:
            raise ValueError("masks_TXY provided but contains no positive pixels")
        center_y, center_x = (ys.min() + ys.max()) // 2, (xs.min() + xs.max()) // 2
        half = min_size_image // 2
        y_start, y_end = center_y - half, center_y - half + min_size_image
        x_start, x_end = center_x - half, center_x - half + min_size_image
        pad_top    = max(0, -y_start)
        pad_bottom = max(0, y_end - image_TZXYC.shape[2])
        pad_left   = max(0, -x_start)
        pad_right  = max(0, x_end - image_TZXYC.shape[3])
        if any((pad_top, pad_bottom, pad_left, pad_right)):
            image_TZXYC = np.pad(image_TZXYC,
                ((0,0),(0,0),(pad_top,pad_bottom),(pad_left,pad_right),(0,0)),
                mode='constant', constant_values=0)
            y_start += pad_top;  y_end += pad_top
            x_start += pad_left; x_end += pad_left

    if cmap_list is None:
        cmap_list = ['gray'] * number_color_channels
    if y_label_list is None:
        y_label_list = [f'Ch {ch}' for ch in range(number_color_channels)]

    fig, ax = plt.subplots(number_color_channels, num_frames,
        figsize=(num_frames * 2, number_color_channels * 2),
        gridspec_kw={'wspace': 0.02, 'hspace': 0.02})
    if number_color_channels == 1:
        ax = np.array([ax])

    for ch in range(number_color_channels):
        for i, frame in enumerate(list_selected_frames):
            current_ax = ax[ch, i]
            if masks_TXY is not None:
                sub = image_TZXYC[frame, 0, y_start:y_end, x_start:x_end, ch]
            else:
                sub = image_TZXYC[frame, 0, :, :, ch]
            current_ax.imshow(sub, vmax=np.percentile(sub, 99.9), cmap=cmap_list[ch])
            if x_title_list is not None and ch == 0:
                current_ax.set_title(x_title_list[i], fontsize=12, fontname='Arial')
            if i == 0:
                current_ax.set_ylabel(y_label_list[ch], fontsize=10, fontname='Arial')
            if coordinates_roi is not None:
                x, y = coordinates_roi[frame]
                if masks_TXY is not None:
                    x -= x_start; y -= y_start
                circ = plt.Circle((x, y), radius_roi_size_px,
                    edgecolor='lightyellow', facecolor='none', linewidth=2)
                current_ax.add_artist(circ)
            sb = ScaleBar(dx=pixel_xy_um, units='um', length_fraction=0.3,
                location='lower right', box_color='black', color='white',
                font_properties={'size': scalebar_size})
            current_ax.add_artist(sb)
            current_ax.set_xticks([]); current_ax.set_yticks([]); current_ax.grid(False)

    if masks_TXY is None and list_axis_limits is not None:
        xmin, xmax, ymin, ymax = list_axis_limits
        for a in ax.flat:
            a.set_xlim(xmin, xmax); a.set_ylim(ymin, ymax)

    plt.tight_layout()
    out_png = results_folder / ('time_courses_' + plot_name + '.png')
    plt.savefig(out_png, dpi=900, bbox_inches='tight', pad_inches=0.1)
    plt.show()


def _compute_kymographs(image_TZXYC, coordinates_roi, length_kymograph_line):
    """Compute kymograph arrays for each channel from ROI coordinates."""
    y = length_kymograph_line
    x = image_TZXYC.shape[0]
    height, width = image_TZXYC.shape[2:4]
    n_channels = image_TZXYC.shape[-1]
    kymographs = []

    for ch in range(n_channels):
        kymograph = np.zeros((y, x))
        for i in range(x):
            cx = np.clip(int(np.round(coordinates_roi[i, 0])), 0, width - 1)
            cy = np.clip(int(np.round(coordinates_roi[i, 1])), 0, height - 1)
            start_x = max(cx - y // 2, 0)
            end_x = min(cx + y // 2, width)
            range_x = np.arange(start_x, end_x, dtype=int)
            line = image_TZXYC[i, 0, cy, range_x, ch]
            kymograph[:len(range_x), i] = line
        kymograph = kymograph[:len(range_x), :]
        kymographs.append(kymograph)
    return kymographs


def plot_kymograph(image_TZXYC, coordinates_roi, list_selected_frames, x_title_list,
                   results_folder, length_kymograph_line=50, cmap_list=None,
                   plot_vertical_lines=False, plot_name='temp'):
    """Generate and plot kymograph from ROI coordinates."""
    kymographs = _compute_kymographs(image_TZXYC, coordinates_roi, length_kymograph_line)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9, 5), sharex=True)
    plt.subplots_adjust(hspace=0, wspace=0)
    for ax, kymo, cmap in zip(axes, kymographs, cmap_list):
        ax.imshow(kymo, aspect='auto', cmap=cmap, vmax=np.percentile(kymo, 99.5))
        ax.grid(False); ax.set_yticks([])
        if plot_vertical_lines:
            for frame in list_selected_frames:
                ax.axvline(x=frame, color='w', linestyle='--', linewidth=4)
        if ax != axes[-1]:
            ax.set_xticks([]); ax.set_xticklabels([])
        else:
            ax.tick_params(axis='x', labelsize=12)
            ax.set_xticks(list_selected_frames)
            ax.set_xticklabels(x_title_list)

    plt.tight_layout()
    plt.savefig(results_folder / f'kymograph_{plot_name}.png', dpi=600)
    plt.savefig(results_folder / f'kymograph_{plot_name}.svg', dpi=600)
    plt.show()


def plot_kymograph_downsampled(image_TZXYC, coordinates_roi, list_selected_frames,
                                x_title_list, list_selected_frame_values_real_time,
                                frame_values, results_folder,
                                length_kymograph_line=50, cmap_list=None,
                                plot_vertical_lines=False, plot_name='temp'):
    """Generate kymograph with temporal downsampling (pre-bleach frames)."""
    # Downsample pre-bleach (1s interval) to match post-bleach (5s interval)
    image_TZXYC_after_pb = image_TZXYC[40:, ...]
    image_TZXYC_before_pb = image_TZXYC[:40:5, ...]
    image_TZXYC_ds = np.concatenate([image_TZXYC_before_pb, image_TZXYC_after_pb], axis=0)

    fv_after = frame_values[39:, ...]
    fv_before = frame_values[:40:5]
    fv_ds = np.concatenate([fv_before, fv_after], axis=0)

    list_selected_frames_ds = [np.argmin(np.abs(fv_ds - x)) for x in list_selected_frame_values_real_time]

    kymographs = _compute_kymographs(image_TZXYC_ds, coordinates_roi, length_kymograph_line)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9, 5), sharex=True)
    plt.subplots_adjust(hspace=0, wspace=0)
    for ax, kymo, cmap in zip(axes, kymographs, cmap_list):
        ax.imshow(kymo, aspect='auto', cmap=cmap, vmax=np.percentile(kymo, 99.5))
        ax.grid(False); ax.set_yticks([])
        if plot_vertical_lines:
            for frame in list_selected_frames_ds:
                ax.axvline(x=frame, color='w', linestyle='--', linewidth=4)
        if ax != axes[-1]:
            ax.set_xticks([]); ax.set_xticklabels([])
        else:
            ax.tick_params(axis='x', labelsize=20)
            ax.set_xticks(list_selected_frames_ds)
            ax.set_xticklabels(x_title_list)

    plt.tight_layout()
    out = results_folder / f'kymograph_ds_{plot_name}.png'
    plt.savefig(out, dpi=900)
    plt.show()
    return out


def plot_merged_image(image_TZXYC, coordinates_roi, length_kymograph_line,
                      cmap_list_imagej, results_folder, pixel_xy_um,
                      list_axis_limits=None, channel_order=[0,1,2],
                      normalize_each_color_channel=False, plot_name='temp',
                      masks_TXY=None, min_size_image=150, gamma=0.6, clip_limit=0.001):
    """Draw merged RGB of one frame with kymograph line, gamma + CLAHE enhancement."""
    frame = 10
    H0, W0 = image_TZXYC.shape[2], image_TZXYC.shape[3]

    if masks_TXY is not None:
        mask_xy = masks_TXY.max(axis=0)
        ys, xs = np.where(mask_xy > 0)
        if ys.size == 0:
            raise ValueError("masks_TXY has no foreground")
        cy, cx = (ys.min()+ys.max())//2, (xs.min()+xs.max())//2
        half = min_size_image // 2
        y0, y1 = cy-half, cy-half+min_size_image
        x0, x1 = cx-half, cx-half+min_size_image
        pad_t, pad_b = max(0,-y0), max(0,y1-H0)
        pad_l, pad_r = max(0,-x0), max(0,x1-W0)
        if any((pad_t, pad_b, pad_l, pad_r)):
            image_TZXYC = np.pad(image_TZXYC,
                ((0,0),(0,0),(pad_t,pad_b),(pad_l,pad_r),(0,0)), mode='constant')
            masks_TXY = np.pad(masks_TXY,
                ((0,0),(pad_t,pad_b),(pad_l,pad_r)), mode='constant')
            y0+=pad_t; y1+=pad_t; x0+=pad_l; x1+=pad_l
        shift = 15
        H_padded = image_TZXYC.shape[2]
        y0 = max(0, y0-shift); y1 = min(H_padded, y1-shift)
        img_crop = image_TZXYC[:, :, y0:y1, x0:x1, :]
        H, W = y1 - y0, x1 - x0
    else:
        img_crop = image_TZXYC
        H, W = H0, W0
        y0 = x0 = 0

    cx = np.clip(int(round(coordinates_roi[frame,0])) - x0, 0, W-1)
    cy = np.clip(int(round(coordinates_roi[frame,1])) - y0, 0, H-1)
    sx, ex = max(cx-length_kymograph_line//2, 0), min(cx+length_kymograph_line//2, W-1)
    range_x = np.arange(sx, ex+1)

    img = img_crop[frame, 0, :, :, :].astype(np.float32)
    C = img.shape[-1]
    merged = np.zeros((H, W, 3), dtype=np.float32)
    for i in range(C):
        ch = img[..., i]
        p1, p99 = np.percentile(ch, 0.1), np.percentile(ch, 99.9)
        ch = np.clip(ch, p1, p99)
        ch = (ch - p1) / (p99 - p1 + 1e-8)
        if normalize_each_color_channel:
            ch = (ch - ch.min()) / (ch.max() - ch.min() + 1e-8)
        ch = gaussian_filter(ch, sigma=0.3)
        ch = ch ** gamma
        col = cmap_list_imagej[i](ch)[..., :3]
        merged += col
    merged /= float(C)
    merged = np.clip(merged, 0, 1)
    merged = exposure.equalize_adapthist(merged, clip_limit=clip_limit)
    if channel_order is not None:
        merged = merged[..., channel_order]

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(merged, interpolation='bicubic')
    ax.plot(range_x, [cy]*len(range_x), color='orangered', linewidth=6, alpha=0.9)
    if masks_TXY is None and list_axis_limits is not None:
        xmin, xmax, ymin, ymax = list_axis_limits
        ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    sb = ScaleBar(dx=pixel_xy_um, units='um', length_fraction=0.3,
                  location='lower right', box_color='black', color='white',
                  font_properties={'size': 0, 'weight': 'bold'})
    ax.add_artist(sb)
    ax.axis('off')
    plt.tight_layout()
    out_png = results_folder / f"{plot_name}.png"
    out_svg = results_folder / f"{plot_name}.svg"
    plt.savefig(out_png, dpi=1500, bbox_inches='tight', pad_inches=0.05,
                facecolor='white', edgecolor='none', format='png')
    plt.savefig(out_svg, bbox_inches='tight', pad_inches=0.05,
                facecolor='white', edgecolor='none', format='svg')
    plt.show()
    return out_png


def save_video_as_avi(image_TZXYC, avi_name, frame_values, pixel_xy_um=0.2,
                      list_axis_limits=None, cmap_list=None, y_label_list=None,
                      fps=5, scalebar_size=0):
    """Export multi-channel microscopy data as an AVI video."""
    n_channels = image_TZXYC.shape[-1]
    if cmap_list is None:
        cmap_list = ['gray'] * n_channels
    ch0_label = y_label_list[0] if y_label_list else None
    ch1_label = y_label_list[1] if y_label_list and len(y_label_list) > 1 else None

    frame_images = []
    for idx, fv in enumerate(frame_values):
        fig, axes = plt.subplots(1, n_channels, figsize=(n_channels*2, 3))
        if n_channels == 1:
            axes = [axes]
        for i in range(n_channels):
            img = image_TZXYC[idx, 0, :, :, i]
            if list_axis_limits is not None:
                x_min, x_max, y_min, y_max = list_axis_limits
                img = img[y_min:y_max, x_min:x_max]
            vmax = np.percentile(img, 99.98) if img.size > 0 else 1
            axes[i].imshow(img, cmap=cmap_list[i], vmax=vmax)
            axes[i].axis('off')
            if i == 0:
                axes[i].text(5, img.shape[0]-5, f"{int(fv)} s",
                    color='white', fontsize=12, ha='left', va='bottom',
                    bbox=dict(facecolor='black', alpha=0.5, pad=2))
        if ch0_label:
            axes[0].text(img.shape[1]//2, 5, ch0_label, color='white', fontsize=10,
                ha='center', va='top', bbox=dict(facecolor='black', alpha=0.5, pad=2))
        if ch1_label:
            axes[1].text(img.shape[1]//2, 5, ch1_label, color='white', fontsize=10,
                ha='center', va='top', bbox=dict(facecolor='black', alpha=0.5, pad=2))
        sb = ScaleBar(dx=pixel_xy_um, units='um', length_fraction=0.3,
            location='lower right', box_color='black', color='white',
            font_properties={'size': scalebar_size})
        axes[-1].add_artist(sb)
        plt.subplots_adjust(wspace=0.05, left=0.01, right=0.99)
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        frame_img = np.array(Image.open(buf))
        frame_img = cv2.cvtColor(frame_img, cv2.COLOR_RGB2BGR)
        frame_images.append(frame_img)

    h, w, _ = frame_images[0].shape
    fourcc = cv2.VideoWriter_fourcc(*'MJPG')
    writer = cv2.VideoWriter(str(avi_name), fourcc, fps, (w, h))
    for frame in frame_images:
        writer.write(frame)
    writer.release()
    print(f"Video saved as {avi_name}")


def compose_pngs(png_path1, png_path2, output_png, spacing=5, target_height=None,
                 bg_color=(255,255,255)):
    """Combine two PNG images side-by-side with high-quality resampling."""
    img1 = Image.open(png_path1).convert("RGB")
    img2 = Image.open(png_path2).convert("RGB")
    if target_height is None:
        target_height = max(img1.height, img2.height)
    target_height = max(target_height, 800)

    def resize(img, h):
        ratio = h / img.height
        new_w = int(img.width * ratio)
        return img.resize((new_w, h), Image.LANCZOS if ratio < 1 else Image.BICUBIC)

    r1, r2 = resize(img1, target_height), resize(img2, target_height)
    canvas = Image.new("RGB", (r1.width + r2.width + spacing, target_height), bg_color)
    canvas.paste(r1, (0, 0))
    canvas.paste(r2, (r1.width + spacing, 0))
    canvas.save(output_png, quality=95, optimize=True)
    canvas.save(output_png.parent / (output_png.stem + "_high_res.png"), quality=100)
    return canvas


# ============================================================================
# Interpretation / Statistics Functions
# ============================================================================

def downsample_to_5_seconds(df):
    """Downsample FRAP data to 5-second intervals (pre-bleach only)."""
    downsampled = []
    for cell_id in df['cell_id'].unique():
        cell_data = df[df['cell_id'] == cell_id].copy()
        high_res = cell_data[cell_data['frame'] < 39]
        low_res = cell_data[cell_data['frame'] >= 39]
        if not high_res.empty:
            downsampled.append(high_res[high_res['frame'] % 5 == 0])
        if not low_res.empty:
            downsampled.append(low_res)
    return pd.concat(downsampled, ignore_index=True) if downsampled else pd.DataFrame()


def plot_FRAP_trajectories(df_list, selected_dataset, results_folder,
                           apply_min_max_normalization=True, display_cell_count=True,
                           bleach_time=10, raw_data_folder=None, figure_panel=None,
                           source_data_condition=None):
    """Plot individual FRAP recovery traces with mean trajectory and one-phase fit.

    Per-cell half-times are fitted only to the post-bleach recovery interval,
    with the bleach time shifted to t=0. When requested, export the individual
    plotted trajectories and plotted time-point means.

    Returns:
        (total_number_cells, t_half_list): Cell count and per-cell t½ values.
    """
    if raw_data_folder is not None and (
        figure_panel is None or source_data_condition is None
    ):
        raise ValueError(
            "figure_panel and source_data_condition are required for source-data export."
        )
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 14, 'axes.labelweight': 'normal'})
    if isinstance(df_list, pd.DataFrame):
        df_list = [df_list]

    def one_phase(t, k):
        return 1 - np.exp(-k * t)

    fig, ax = plt.subplots(figsize=(6, 4), facecolor='white')
    ax.set_facecolor('white')
    all_data, source_data, t_half_list = [], [], []
    total_cells = 0
    col = 'mean_roi_frap'

    for df in df_list:
        df_sel = df[df['dataset_type'] == selected_dataset]
        for _, row in df_sel[['subfolder_id', 'image_name']].drop_duplicates().iterrows():
            cell_data = df_sel[
                (df_sel['subfolder_id'] == row['subfolder_id']) &
                (df_sel['image_name'] == row['image_name'])
            ].copy()
            if apply_min_max_normalization:
                mn, mx = cell_data[col].min(), cell_data[col].max()
                cell_data[col] = (cell_data[col] - mn) / (mx - mn) if mx > mn else 0.0
            try:
                recovery_data = cell_data[cell_data['frame'] >= bleach_time]
                recovery_time = recovery_data['frame'].values - bleach_time
                popt, _ = curve_fit(one_phase, recovery_time,
                    recovery_data[col].values, p0=[0.1], bounds=(0, np.inf))
                t_half_list.append(np.log(2) / popt[0])
            except Exception:
                t_half_list.append(np.nan)
            ax.plot(cell_data['frame'], cell_data[col], '-', color='dimgray', lw=0.5, alpha=0.2)
            all_data.append(cell_data[['frame', col]])
            total_cells += 1
            if raw_data_folder is not None:
                cell_source_data = cell_data[['frame', col]].rename(columns={
                    'frame': 'time_s',
                    col: 'normalized_mean_roi_intensity',
                })
                cell_source_data.insert(0, 'cell_id', f'cell_{total_cells:03d}')
                cell_source_data.insert(0, 'condition', source_data_condition)
                source_data.append(cell_source_data)

    if all_data:
        all_df = pd.concat(all_data, ignore_index=True)
        summary_data = all_df.groupby('frame')[col].agg(['mean', 'count']).reset_index()
        ax.plot(
            summary_data['frame'], summary_data['mean'], '-', color='green',
            lw=3, label='Mean Trajectory',
        )
        if raw_data_folder is not None:
            raw_data_folder.mkdir(parents=True, exist_ok=True)
            source_data_df = pd.concat(source_data, ignore_index=True)
            source_data_path = (
                raw_data_folder /
                f'Figure_{figure_panel}_{source_data_condition}_FRAP_individual_trajectories.csv'
            )
            source_data_df.to_csv(source_data_path, index=False)

            summary_data = summary_data.rename(columns={
                'frame': 'time_s',
                'mean': 'mean_normalized_intensity',
                'count': 'n_cells',
            })
            summary_data.insert(0, 'condition', source_data_condition)
            summary_data_path = (
                raw_data_folder /
                f'Figure_{figure_panel}_{source_data_condition}_FRAP_timepoint_summary.csv'
            )
            summary_data.to_csv(summary_data_path, index=False)

    ax.set_xlabel("Time (s)", fontdict={'family': 'Arial', 'size': 20, 'color': 'black'})
    ax.set_ylabel("Normalized Intensity", fontdict={'family': 'Arial', 'size': 20, 'color': 'black'})
    ax.tick_params(axis='both', which='major', labelsize=16, colors='black')
    for spine in ax.spines.values():
        spine.set_color('black')
    ax.grid(False)
    if display_cell_count:
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.text(xlim[1]*0.95, ylim[0]+0.05*(ylim[1]-ylim[0]),
            f"n = {total_cells}", ha='right', va='bottom', fontsize=20, fontname="Arial", color='black')
    fig.tight_layout()
    fname = f"{selected_dataset}_FRAP_trajectories.png"
    if selected_dataset == 'UTag($\\Delta$Cys)':
        fname = "UTag_deltaC_FRAP_trajectories.png"
    fpath = results_folder / fname
    fig.savefig(fpath, dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    fig.savefig(fpath.with_suffix('.svg'), format='svg', bbox_inches='tight', pad_inches=0.1, transparent=True)
    plt.show()
    return total_cells, t_half_list


def save_t_half_source_data(t_half_lists, conditions, raw_data_folder):
    """Save the per-cell fitted half-times plotted in Figure 4H."""
    if len(t_half_lists) != len(conditions):
        raise ValueError("t_half_lists and conditions must have the same length.")
    source_data = []
    for condition, t_half_values in zip(conditions, t_half_lists):
        condition_data = pd.DataFrame({
            'condition': condition,
            't_half_s': np.asarray(t_half_values),
        }).dropna(subset=['t_half_s'])
        source_data.append(condition_data)
    source_data_df = pd.concat(source_data, ignore_index=True)
    raw_data_folder.mkdir(parents=True, exist_ok=True)
    source_data_path = raw_data_folder / 'Figure_4H_FRAP_t_half_individual_values.csv'
    source_data_df.to_csv(source_data_path, index=False)
    return source_data_df


def write_frap_source_data_readme(raw_data_folder):
    """Write the Figure 4 source-data dictionary and analysis definitions."""
    readme_text = """# Figure 4 FRAP source data

These files contain the processed numerical data used to generate the graphs in
Figure 4A-H. They are figure source data, not unprocessed acquisition-level raw
microscopy data. Files are regenerated from the in-memory data assembled from
the experiment-level CSV files in `FRAP_quantification/`;
`combined_FRAP_data.csv` is not read or required.

## Condition labels

| CSV label | Figure label |
|---|---|
| `UTag` | UTag |
| `UTag_DeltaCys` | UTag(ΔCys), ΔCys |
| `SunTag` | SunTag |
| `ALFA_tag` | ALFA-tag, ALFA |
| `HA` | HA |

ASCII-only condition labels are used in every CSV to avoid encoding and naming
inconsistencies.

## Figure 4A-E: condition-specific FRAP trajectories

Each panel has an individual-trajectory file and a time-point summary file.
The individual files contain `condition`, anonymous `cell_id`, `time_s`, and
`normalized_mean_roi_intensity`. Each value is the mean fluorescence intensity
within one cell ROI at one time point after per-cell min-max normalization:
`(value - cell minimum) / (cell maximum - cell minimum)`.

The summary files contain `condition`, `time_s`,
`mean_normalized_intensity`, and `n_cells`. The green line is the arithmetic
mean across the available cells at each time point. Panels A-E display no error
bars or error bands, so no variability column is included.

| Panel | Condition | Individual rows | Cells |
|---|---|---:|---:|
| 4A | `UTag` | 7,000 | 50 |
| 4B | `UTag_DeltaCys` | 6,160 | 44 |
| 4C | `SunTag` | 7,980 | 57 |
| 4D | `ALFA_tag` | 7,980 | 57 |
| 4E | `HA` | 5,880 | 42 |

## Figure 4F: combined FRAP trajectories

`Figure_4F_FRAP_all_conditions_timepoint_summary.csv` contains `condition`,
`time_s`, `mean_normalized_intensity`,
`sample_sd_normalized_intensity`, and `n_cells`. The colored line is the
arithmetic mean across cells and the shaded region is mean ± sample standard
deviation (SD; pandas default `ddof=1`). It is not SEM.

The existing Figure 4F plotting call repeats the recovery-drop quality check
with a threshold of 0.8 on the prepared normalized trajectories. It retains
49 UTag, 44 UTag_DeltaCys, 57 SunTag, 57 ALFA_tag, and 42 HA cells. This is why
the UTag `n_cells` in panel F is 49 rather than the 50 used in panel A.

## Figure 4G: recovered intensity

`Figure_4G_FRAP_recovered_intensity_individual_values.csv` contains
`condition` and `normalized_recovered_intensity`. Each row is one black dot:
the final recorded normalized mean ROI intensity for one cell. Box medians,
quartiles, whiskers, statistical tests, p-values, and significance labels are
not included.

## Figure 4H: recovery half-time

`Figure_4H_FRAP_t_half_individual_values.csv` contains `condition` and
`t_half_s`. Each row is one successfully fitted black dot. Fits use only the
post-bleach interval, with the 10-second bleach point reset to t=0. The model is
`I(t) = 1 - exp(-k*t)`, with the recovery plateau fixed at 1, and
`t_half_s = ln(2) / k`. Box summaries and statistical-test results are not
included.

## Data preparation

The initial notebook quality check retains trajectories whose fluorescence drop
during the first 20 seconds is greater than 0.5. The files contain 250 cells
from 15 experiment folders: 50 UTag from 3 folders, 44 UTag_DeltaCys from 3,
57 SunTag from 4, 57 ALFA_tag from 2, and 42 HA from 3. Figure 4F has the
additional filter described above.

There are 140 time points per complete trajectory: 0-39 seconds in one-second
increments, followed by 44-539 seconds in five-second increments. Optional
downsampling is disabled. The final values in Figure 4G are at 539 seconds.

Normalization divisions are guarded: a constant trajectory is assigned zero
after normalization. The CSV exports are created from the same prepared values
used by their corresponding plotting functions. The package contains no
p-values, significance tests, box summaries, image-provenance fields, source
image names, experiment identifiers, or replacement combined trajectory file.
"""
    raw_data_folder.mkdir(parents=True, exist_ok=True)
    readme_path = raw_data_folder / 'README_source_data.md'
    readme_path.write_text(readme_text, encoding='utf-8')
    return readme_path


def plot_t_half_box_swarm(t_half_lists, names, results_folder, figsize=(6, 4),
                          ylabel="t₁/₂ (s)", title="", y_min=None, y_max=None,
                          swarm_color="black", tick_size=16, show_stats=False):
    """Plot boxplots + swarm of per-trajectory t½ values with optional pairwise statistics."""
    lengths = [len(lst) for lst in t_half_lists]
    df = pd.DataFrame({
        "Condition": np.repeat(names, lengths),
        "t_half": np.concatenate([np.asarray(lst) for lst in t_half_lists])
    })
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    sns.boxplot(x="Condition", y="t_half", data=df, order=names, showfliers=False, ax=ax,
        boxprops={'facecolor': 'white', 'edgecolor': 'black'},
        medianprops={'color': 'red'}, whiskerprops={'color': 'black'}, capprops={'color': 'black'})
    sns.swarmplot(x="Condition", y="t_half", data=df, order=names, color=swarm_color, ax=ax)
    ax.set_xlabel("", fontsize=tick_size, fontname="Arial", color="black")
    ax.set_ylabel(ylabel, fontsize=tick_size+4, fontname="Arial", color="black")
    if title:
        ax.set_title(title, fontsize=tick_size+4, fontname="Arial", color="black")
    if y_min is not None or y_max is not None:
        ax.set_ylim(bottom=y_min, top=y_max)
    ax.tick_params(axis='x', labelsize=tick_size+2, colors='black')
    ax.tick_params(axis='y', labelsize=tick_size-2, colors='black')
    for label in ax.get_xticklabels():
        label.set_fontfamily("Arial")
    for label in ax.get_yticklabels():
        label.set_fontfamily("Arial")

    if show_stats:
        combos = list(combinations(range(len(names)), 2))
        y_lo, y_hi = ax.get_ylim()
        offset = 0.05 * (y_hi - y_lo)
        for level, (i, j) in enumerate(combos, start=1):
            g1 = df.loc[df["Condition"]==names[i], "t_half"]
            g2 = df.loc[df["Condition"]==names[j], "t_half"]
            _, p = stats.mannwhitneyu(g1, g2)
            sig = "****" if p<1e-4 else "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns"
            y_line = max(g1.max(), g2.max()) + offset * level
            h = offset * 0.05
            ax.plot([i,i,j,j], [y_line, y_line+h, y_line+h, y_line], lw=1.5, c='k')
            ax.text((i+j)/2, y_line+h, sig, ha='center', va='bottom',
                fontsize=tick_size, fontname="Arial", color="k")
            print(f"Comparison {names[i]} vs {names[j]}: p = {p:.4f}")

    fig.tight_layout()
    fpath = results_folder / "t_half_box_swarm.png"
    fig.savefig(fpath, dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    fig.savefig(fpath.with_suffix('.svg'), format='svg', bbox_inches='tight', pad_inches=0.1, transparent=True)
    plt.show()
    return ax


def plot_mean_trajectories_all(df, selected_datasets, results_folder,
                                selected_field='mean_roi_frap',
                                apply_quality_check=True, drop_threshold=0.2,
                                fontsize=12, apply_min_max_normalization=True,
                                fig_size=(6,4), use_sem=True, raw_data_folder=None,
                                condition_label_map=None):
    """Plot normalized mean FRAP recovery curves with SD or SEM shading."""
    cmap = plt.get_cmap('tab20')
    color_map = {ds: cmap(i) for i, ds in enumerate(selected_datasets)}

    fig, ax = plt.subplots(figsize=fig_size)
    ax.set_facecolor('white')
    source_data = []
    for ds in selected_datasets:
        df_ds = df[df['dataset_type'] == ds]
        if df_ds.empty:
            continue
        curves = []
        for (sub, img), group in df_ds.groupby(['subfolder_id', 'image_name']):
            series = group[['frame', selected_field]].copy()
            if apply_quality_check:
                init = series[selected_field].iloc[0]
                drop = init - series.loc[series['frame'] <= 20, selected_field].min()
                if drop <= drop_threshold:
                    continue
            if apply_min_max_normalization:
                mn, mx = series[selected_field].min(), series[selected_field].max()
                if mx > mn:
                    series[selected_field] = (series[selected_field] - mn) / (mx - mn)
                else:
                    series[selected_field] = 0
            curves.append(series.set_index('frame'))
        if not curves:
            continue
        all_cells = pd.concat(curves, axis=1)
        means = all_cells.mean(axis=1)
        errs = all_cells.sem(axis=1) if use_sem else all_cells.std(axis=1)
        c = color_map[ds]
        ax.plot(means.index, means.values, 'o-', color=c, lw=2, label=ds, markersize=3)
        ax.fill_between(means.index, means - errs, means + errs, color=c, alpha=0.1)
        if raw_data_folder is not None:
            if use_sem:
                error_column = 'sem_normalized_intensity'
            else:
                error_column = 'sample_sd_normalized_intensity'
            condition = condition_label_map.get(ds, ds) if condition_label_map else ds
            condition_source_data = pd.DataFrame({
                'condition': condition,
                'time_s': means.index,
                'mean_normalized_intensity': means.values,
                error_column: errs.values,
                'n_cells': all_cells.count(axis=1).values,
            })
            source_data.append(condition_source_data)

    if raw_data_folder is not None and source_data:
        raw_data_folder.mkdir(parents=True, exist_ok=True)
        source_data_df = pd.concat(source_data, ignore_index=True)
        source_data_path = (
            raw_data_folder /
            'Figure_4F_FRAP_all_conditions_timepoint_summary.csv'
        )
        source_data_df.to_csv(source_data_path, index=False)

    ax.set_ylim(-0.05, 1.5)
    for spine in ax.spines.values():
        spine.set_color('black'); spine.set_linewidth(1)

    font_prop = FontProperties(family='Arial', size=fontsize)
    leg = ax.legend(fontsize=fontsize-1, loc='upper center', bbox_to_anchor=(0.5, 1.4),
        ncol=len(selected_datasets), frameon=True, framealpha=1, edgecolor='black',
        facecolor='white', prop=font_prop, columnspacing=1.0, handletextpad=1)
    plt.setp(leg.get_texts(), fontfamily='Arial', color='black')
    ax.tick_params(axis='both', which='major', labelsize=fontsize+1, colors='black')
    fig.subplots_adjust(right=0.75)

    out_png = results_folder / "mean_FRAP_trajectories_all.png"
    fig.savefig(out_png, dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    fig.savefig(out_png.with_suffix('.svg'), format='svg', bbox_inches='tight', pad_inches=0.1, transparent=True)
    plt.show()
    return ax


def plot_box_swarm_final_values(df, selected_field, results_folder, figsize=(6, 4),
                                 xlabel="Dataset Type", ylabel="Final Normalized Intensity",
                                 title="", y_min=None, y_max=None, swarm_color="black",
                                 tick_size=16, show_stats=False, raw_data_folder=None,
                                 condition_label_map=None):
    """Boxplot + swarmplot of final-frame FRAP values per cell, with optional pairwise stats."""
    sns.set_style("ticks")
    order_categories = ['UTag', 'UTag($\\Delta$Cys)', 'SunTag', 'ALFA-tag', 'HA']
    final_df = df.loc[df.groupby('cell_id')['frame'].idxmax()].copy()
    if raw_data_folder is not None:
        source_data_df = final_df[['dataset_type', selected_field]].rename(columns={
            'dataset_type': 'condition',
            selected_field: 'normalized_recovered_intensity',
        })
        if condition_label_map:
            source_data_df['condition'] = source_data_df['condition'].replace(
                condition_label_map
            )
        source_data_df = source_data_df[
            ['condition', 'normalized_recovered_intensity']
        ]
        raw_data_folder.mkdir(parents=True, exist_ok=True)
        source_data_path = (
            raw_data_folder /
            'Figure_4G_FRAP_recovered_intensity_individual_values.csv'
        )
        source_data_df.to_csv(source_data_path, index=False)

    fig, ax = plt.subplots(figsize=figsize, facecolor='white')
    sns.boxplot(x="dataset_type", y=selected_field, data=final_df,
        order=order_categories, showfliers=False, ax=ax,
        boxprops={'facecolor': 'white', 'edgecolor': 'black'},
        medianprops={'color': 'red'}, whiskerprops={'color': 'black'}, capprops={'color': 'black'})
    ax.set_facecolor('white')
    sns.swarmplot(x="dataset_type", y=selected_field, data=final_df,
        order=order_categories, color=swarm_color, ax=ax)

    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=tick_size+4, fontname="Arial", color='black')
    if y_min is not None and y_max is not None:
        ax.set_ylim(y_min, y_max)
    ax.tick_params(axis='x', labelsize=tick_size+4, colors='black')
    ax.tick_params(axis='y', labelsize=tick_size, colors='black')
    for label in ax.get_xticklabels():
        label.set_fontfamily("Arial")
    for label in ax.get_yticklabels():
        label.set_fontfamily("Arial")

    if show_stats:
        comparisons = [((i, j), j - i) for i, j in combinations(range(len(order_categories)), 2)]
        y_lo, y_hi = ax.get_ylim()
        offset = 0.1 * (y_hi - y_lo)
        for ((i, j), level) in comparisons:
            g1 = final_df.loc[final_df['dataset_type']==order_categories[i], selected_field]
            g2 = final_df.loc[final_df['dataset_type']==order_categories[j], selected_field]
            _, p = stats.mannwhitneyu(g1, g2)
            sig = '****' if p<1e-4 else '***' if p<1e-3 else '**' if p<1e-2 else '*' if p<0.05 else 'ns'
            y_line = max(g1.max(), g2.max()) + offset * level
            h = offset * 0.2
            ax.plot([i,i,j,j], [y_line, y_line+h, y_line+h, y_line], lw=1.5, c='k')
            ax.text((i+j)*.5, y_line+h, sig, ha='center', va='bottom',
                fontsize=tick_size, fontname="Arial", color='k')
            print(f"Comparison {order_categories[i]} vs {order_categories[j]}: p = {p:.6f}")

    ax.set_xticklabels(order_categories, fontsize=tick_size, fontname="Arial", color='black')
    fig.tight_layout()

    # Report statistics
    for cat in order_categories:
        subset = final_df[final_df['dataset_type'] == cat][selected_field]
        mean_v, std_v = subset.mean(), subset.std()
        sem_v = std_v / np.sqrt(len(subset)) if len(subset) > 0 else 0
        print(f"{cat}: mean = {mean_v:.2f}, std = {std_v:.2f}, sem = {sem_v:.2f}")

    fpath = results_folder / f"box_swarm_{selected_field}.png"
    fig.savefig(fpath, dpi=300, bbox_inches='tight', pad_inches=0.1, transparent=True)
    fig.savefig(fpath.with_suffix('.svg'), format='svg', bbox_inches='tight', pad_inches=0.1, transparent=True)
    plt.show()
    return ax
