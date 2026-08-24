"""Self-contained synthetic translation spot size vs. intensity comparison.

Purpose
-------
Illustration for the reviewer question on why SunTag translation sites can appear
LARGER (greater Gaussian FWHM) without showing a correspondingly HIGHER measured
spot intensity or signal-to-noise ratio.

Key idea
--------
"Spot size" and "spot intensity" are *independent* measurements in the pipeline:

  * Spot size      = FWHM of the fitted 2D Gaussian  = 2*sqrt(2*ln2) * sigma
  * Spot intensity = mean(disk) - mean(doughnut)      (disk-doughnut, LOCAL background)
  * SNR (mean)     = (mean(disk) - mean(doughnut)) / std(doughnut)

We construct two spots with DIFFERENT sigmas (hence different FWHM / size) that are
forced to share BOTH an identical measured disk-doughnut intensity AND an identical
mean-based SNR, while the local background level is left free to differ. This needs
two degrees of freedom:

  1. The wide-spot peak amplitude is solved so the two spots have the same measured
     disk-doughnut intensity.
  2. Each spot is given its own background-noise level (its "free" local background),
     solved so the two spots have the same doughnut std and therefore the same SNR.

The result directly mirrors Figures 6C-D: comparable spot intensity and SNR across
reporters, with a larger apparent spot size for the broader (SunTag-like) reporter.

Faithfulness to the analysis pipeline (microlive/microlive/microscopy.py)
-------------------------------------------------------------------------
The disk / doughnut geometry and the intensity / SNR formulas below replicate the
SpotIntensity.calculate_intensity() implementation:

  * disk      : square of side  best_size = spot_size + 2   (microscopy.py:2221, 2333-2343)
  * doughnut  : bounded ring -> crop of side best_size + PIXELS_AROUND_SPOT, with the
                central best_size x best_size square removed   (microscopy.py:2007, 2148-2156, 2344-2349)
  * intensity : mean(disk) - mean(doughnut)                   (microscopy.py:2185-2190)
  * SNR       : 'peak'          -> (max(disk)  - mean(donut)) / std(donut)
                'disk_doughnut' -> (mean(disk) - mean(donut)) / std(donut)  (microscopy.py:2173-2182)
  * spot size : 2*sqrt(2*ln2) * sigma  (FWHM)                 (microscopy.py:6029-6037)

The background here is the MEASURED doughnut mean (not a fixed baseline), so the
width-dependent leakage of the Gaussian tail into the doughnut is reproduced.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Patch
from matplotlib.colors import ListedColormap
from pathlib import Path
from scipy.optimize import brentq

# ---------------------------------------------------------------------
# Global text style: Arial, black, non-bold, consistent sizes everywhere
# ---------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
    "font.weight": "normal",
    "axes.titleweight": "normal",
    "axes.labelweight": "normal",
    "text.color": "black",
    "axes.labelcolor": "black",
    "axes.titlecolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
})
TITLE_SIZE = 12   # panel titles
LABEL_SIZE = 10   # axis labels (X / Y / Z / Intensity)
TICK_SIZE = 8     # tick numbers
ANNOT_SIZE = 10   # in-panel annotations and legend

# ---------------------------------------------------------------------
# User-relevant parameters
# ---------------------------------------------------------------------
CROP_SIZE = 13             # synthetic field of view (px); large enough for the doughnut crop
SIGMA_NARROW = 1.3        # Gaussian sigma of the narrow spot (typical / UTag-like)
SIGMA_WIDE = 1.5          # Gaussian sigma of the wide spot (longer array / SunTag-like)
BASELINE = 100.0          # true image background (used only to build the synthetic image)
AMPLITUDE_NARROW = 300.0  # peak amplitude of the narrow spot (sets the brightness scale)
SNR_TARGET = 5.0          # common mean-based SNR enforced on BOTH spots
SEED = 7
SPOT_SIZE = 5             # pipeline `spot_size` (default 5)
PIXELS_AROUND_SPOT = 3    # pipeline PIXELS_AROUND_SPOT (microscopy.py:2007)
SNR_METHOD = "disk_doughnut"  # 'disk_doughnut' = mean-based SNR; 'peak' uses max (microscopy.py:2173-2182)

# Pipeline derives the disk side as spot_size + 2 ("+2 to ensure full spot coverage")
best_size = SPOT_SIZE + 2

# Output configuration
results_folder = Path(__file__).parent / "results_spots_synthetic"
results_folder.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------
center = (CROP_SIZE - 1) / 2
center_idx = int(round(center))
x = np.arange(CROP_SIZE)
y = np.arange(CROP_SIZE)
X, Y = np.meshgrid(x, y)
R2 = (X - center) ** 2 + (Y - center) ** 2


# ---------------------------------------------------------------------
# Disk / doughnut intensity measurement (faithful to microscopy.py)
# ---------------------------------------------------------------------
def return_crop(image, cx, cy, spot_range):
    """Replicates SpotIntensity.calculate_intensity.return_crop (microscopy.py:2143-2146)."""
    cx, cy = int(cx), int(cy)
    return image[cy + spot_range[0]:cy + (spot_range[-1] + 1),
                 cx + spot_range[0]:cx + (spot_range[-1] + 1)].copy()


def return_donut(crop, spot_size):
    """Replicates SpotIntensity.calculate_intensity.return_donut (microscopy.py:2148-2156).

    Removes the central spot_size x spot_size square; the remaining bounded ring is
    the local background. (The pipeline casts to uint16; we keep float for clarity,
    which does not change the conceptual result.)
    """
    tem = crop.copy().astype(float)
    c = tem.shape[0] // 2
    half = spot_size // 2
    tem[c - half:c + half + 1, c - half:c + half + 1] = np.nan
    return tem[~np.isnan(tem)]


def measure_disk_doughnut(img, cx, cy, snr_method=SNR_METHOD):
    """Measure spot intensity and SNR using the disk-doughnut pipeline geometry."""
    half = best_size // 2
    values_disk = img[cy - half:cy + half + 1, cx - half:cx + half + 1]

    crop_size = best_size + PIXELS_AROUND_SPOT
    crop_range = np.arange(-(crop_size - 1) / 2, (crop_size - 1) / 2 + 1, 1, dtype=int)
    crop_disk_and_donut = return_crop(img, cx, cy, crop_range)
    values_donut = return_donut(crop_disk_and_donut, spot_size=best_size)

    mean_disk = values_disk.mean()
    mean_donut = values_donut.mean()
    std_donut = values_donut.std()

    intensity = mean_disk - mean_donut                          # disk_donut (microscopy.py:2189)
    snr_peak = (values_disk.max() - mean_donut) / std_donut if std_donut > 0 else 0.0
    snr_disk_donut = intensity / std_donut if std_donut > 0 else 0.0
    snr = snr_peak if snr_method == "peak" else snr_disk_donut

    return {
        "intensity": intensity,
        "snr": snr,
        "snr_peak": snr_peak,
        "snr_disk_donut": snr_disk_donut,
        "mean_disk": mean_disk,
        "mean_donut": mean_donut,
        "std_donut": std_donut,
        "peak": values_disk.max(),
        "total": values_disk.sum() - mean_donut * values_disk.size,  # background-subtracted integrated signal
    }


def clean_intensity_coeff(sigma):
    """Clean (noise-free) disk-doughnut intensity for a unit-amplitude Gaussian.

    intensity(A, sigma) = A * coeff(sigma), because the constant baseline cancels in
    (mean_disk - mean_donut).
    """
    unit_spot = np.exp(-R2 / (2 * sigma ** 2))  # A=1, baseline=0, no noise
    return measure_disk_doughnut(unit_spot, center_idx, center_idx)["intensity"]


# ---------------------------------------------------------------------
# Disk / doughnut pixel masks (exactly the pixels measure_disk_doughnut reads)
# ---------------------------------------------------------------------
def build_masks():
    half = best_size // 2
    disk_mask = np.zeros((CROP_SIZE, CROP_SIZE), dtype=bool)
    disk_mask[center_idx - half:center_idx + half + 1,
              center_idx - half:center_idx + half + 1] = True

    crop_size = best_size + PIXELS_AROUND_SPOT
    crop_range = np.arange(-(crop_size - 1) / 2, (crop_size - 1) / 2 + 1, 1, dtype=int)
    r0 = center_idx + crop_range[0]
    r1 = center_idx + crop_range[-1] + 1
    cw = r1 - r0
    c = cw // 2
    hole0 = r0 + (c - half)
    hole1 = r0 + (c + half + 1)
    crop_mask = np.zeros((CROP_SIZE, CROP_SIZE), dtype=bool)
    crop_mask[r0:r1, r0:r1] = True
    hole_mask = np.zeros((CROP_SIZE, CROP_SIZE), dtype=bool)
    hole_mask[hole0:hole1, hole0:hole1] = True
    donut_mask = crop_mask & ~hole_mask
    return disk_mask, donut_mask


DISK_MASK, DONUT_MASK = build_masks()


def demean_regions(z, masks, iters=500, tol=1e-13):
    """Remove the mean of z over each mask (alternating projection for overlaps).

    After this, mean(z[disk]) = mean(z[donut]) = 0, so scaled noise contributes
    exactly zero to the measured disk-doughnut intensity -> the intensity equals the
    clean (noise-free) value regardless of the noise amplitude.
    """
    z = z.copy()
    for _ in range(iters):
        worst = 0.0
        for m in masks:
            mu = z[m].mean()
            z[m] -= mu
            worst = max(worst, abs(mu))
        if worst < tol:
            break
    return z


def solve_noise_sd(clean_signal, z, target_std):
    """Solve the noise amplitude so the measured doughnut std equals target_std."""
    def f(sd):
        img = BASELINE + clean_signal + sd * z
        return measure_disk_doughnut(img, center_idx, center_idx)["std_donut"] - target_std
    # std at sd=0 is the signal-gradient floor; target must exceed it (adding noise only raises std)
    if f(0.0) > 0:
        raise ValueError("SNR_TARGET too high for this sigma: signal-gradient noise alone exceeds target.")
    return brentq(f, 0.0, 1000.0)


# ---------------------------------------------------------------------
# Construction: same intensity AND same SNR, different sigma/size, free background
# ---------------------------------------------------------------------
# Knob 1 -- match the measured disk-doughnut intensity:
clean_intensity = AMPLITUDE_NARROW * clean_intensity_coeff(SIGMA_NARROW)
AMPLITUDE_WIDE = clean_intensity / clean_intensity_coeff(SIGMA_WIDE)

# Knob 2 -- match the SNR: choose a common doughnut std, then solve each spot's
# background-noise level (its free local background) to reach it.
target_std = clean_intensity / SNR_TARGET

clean_narrow = AMPLITUDE_NARROW * np.exp(-R2 / (2 * SIGMA_NARROW ** 2))
clean_wide = AMPLITUDE_WIDE * np.exp(-R2 / (2 * SIGMA_WIDE ** 2))

rng = np.random.default_rng(SEED)
z_narrow = demean_regions(rng.standard_normal((CROP_SIZE, CROP_SIZE)), (DISK_MASK, DONUT_MASK))
z_wide = demean_regions(rng.standard_normal((CROP_SIZE, CROP_SIZE)), (DISK_MASK, DONUT_MASK))

noise_sd_narrow = solve_noise_sd(clean_narrow, z_narrow, target_std)
noise_sd_wide = solve_noise_sd(clean_wide, z_wide, target_std)

image_narrow = BASELINE + clean_narrow + noise_sd_narrow * z_narrow
image_wide = BASELINE + clean_wide + noise_sd_wide * z_wide

m_narrow = measure_disk_doughnut(image_narrow, center_idx, center_idx)
m_wide = measure_disk_doughnut(image_wide, center_idx, center_idx)

# ---------------------------------------------------------------------
# Gaussian FWHM (spot size) -- microscopy.py:6029-6037
# ---------------------------------------------------------------------
fwhm_factor = 2 * np.sqrt(2 * np.log(2))
fwhm_narrow = fwhm_factor * SIGMA_NARROW
fwhm_wide = fwhm_factor * SIGMA_WIDE

# ---------------------------------------------------------------------
# Figure: 2x2 -- top row 2D disk/doughnut grids; bottom row 3D surfaces
# ---------------------------------------------------------------------
max_z = max(image_narrow.max(), image_wide.max()) + 50

# Illustrative disk (signal) and doughnut (background) regions, drawn concentric
# for clarity. The disk (7x7) is exactly the pipeline's; the doughnut ring is shown
# concentric -- the actual crop is marginally asymmetric (microscopy.py:2344-2349).
half_disk = best_size // 2
half_outer = (best_size + PIXELS_AROUND_SPOT) // 2
ILL_DISK = (np.abs(X - center) <= half_disk) & (np.abs(Y - center) <= half_disk)
ILL_DONUT = ((np.abs(X - center) <= half_outer) & (np.abs(Y - center) <= half_outer)) & ~ILL_DISK

DISK_COLOR = (0.85, 0.12, 0.12)
DONUT_COLOR = (0.10, 0.55, 0.50)


def plot_disk_doughnut_grid(ax, image, sigma, label):
    """2D top view of the crop with the 13x13 pixel grid and disk/doughnut regions."""
    ext = [-0.5, CROP_SIZE - 0.5, -0.5, CROP_SIZE - 0.5]
    ax.imshow(image, cmap="coolwarm", origin="lower", vmin=0, vmax=max_z, extent=ext, zorder=1)

    ax.add_patch(Rectangle((center_idx - half_disk - 0.5, center_idx - half_disk - 0.5),
                           2 * half_disk + 1, 2 * half_disk + 1,
                           fill=False, edgecolor=DISK_COLOR, lw=2, zorder=3))
    ax.add_patch(Rectangle((center_idx - half_outer - 0.5, center_idx - half_outer - 0.5),
                           2 * half_outer + 1, 2 * half_outer + 1,
                           fill=False, edgecolor=DONUT_COLOR, lw=2, zorder=3))

    # 13x13 pixel grid: labelled ticks at integer pixels, gridlines on cell boundaries
    ax.set_xticks(np.arange(0, CROP_SIZE))
    ax.set_yticks(np.arange(0, CROP_SIZE))
    ax.set_xticks(np.arange(-0.5, CROP_SIZE, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, CROP_SIZE, 1), minor=True)
    ax.grid(which="minor", color="0.55", lw=0.6, zorder=4)
    ax.tick_params(labelsize=TICK_SIZE, colors="black")
    ax.set_xlim(-0.5, CROP_SIZE - 0.5)
    ax.set_ylim(-0.5, CROP_SIZE - 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("X (px)", fontsize=LABEL_SIZE, color="black")
    ax.set_ylabel("Y (px)", fontsize=LABEL_SIZE, color="black")
    ax.set_title(rf"{label} ($\sigma={sigma:.1f}$)", fontsize=TITLE_SIZE, color="black")


def style_axis(ax, image, measurement, fwhm):
    surf = ax.plot_surface(X, Y, image, cmap="coolwarm",
                           edgecolor="none", vmin=0, vmax=max_z)
    # Intensity and SNR are matched between panels; size and background differ.
    annotation = (f"Intensity: {measurement['intensity']:.1f}\n"
                  f"SNR: {measurement['snr']:.1f}\n"
                  f"Spot size: {fwhm:.2f} px\n"
                  f"Background: {measurement['mean_donut']:.1f}")
    ax.text2D(0.05, 0.80, annotation, transform=ax.transAxes, fontsize=ANNOT_SIZE, color="black")
    ax.set_xlim(0, CROP_SIZE - 1)
    ax.set_ylim(0, CROP_SIZE - 1)
    ax.set_zlim(0, max_z)
    ax.set_box_aspect((1, 1, 0.8))
    ax.set_xlabel("X (px)", fontsize=LABEL_SIZE, color="black")
    ax.set_ylabel("Y (px)", fontsize=LABEL_SIZE, color="black")
    ax.set_zlabel("Intensity", fontsize=LABEL_SIZE, color="black")
    ax.tick_params(axis="both", labelsize=TICK_SIZE, colors="black")
    ax.tick_params(axis="z", labelsize=TICK_SIZE, colors="black")
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.xaxis._axinfo["grid"]["color"] = (0.7, 0.7, 0.7, 0.15)
    ax.yaxis._axinfo["grid"]["color"] = (0.7, 0.7, 0.7, 0.15)
    ax.zaxis._axinfo["grid"]["color"] = (0.7, 0.7, 0.7, 0.15)
    ax.view_init(elev=25, azim=-45)
    return surf


fig = plt.figure(figsize=(10, 10))
gs = GridSpec(2, 2, figure=fig, height_ratios=[0.85, 1.1],
              hspace=0.26, wspace=0.22, left=0.08, right=0.95, top=0.93, bottom=0.12)

# Top row: 2D disk/doughnut pixel grids
ax_2d_n = fig.add_subplot(gs[0, 0])
ax_2d_w = fig.add_subplot(gs[0, 1])
plot_disk_doughnut_grid(ax_2d_n, image_narrow, SIGMA_NARROW, "Narrow")
plot_disk_doughnut_grid(ax_2d_w, image_wide, SIGMA_WIDE, "Wide")

# Bottom row: 3D surfaces
ax_3d_n = fig.add_subplot(gs[1, 0], projection="3d")
surf1 = style_axis(ax_3d_n, image_narrow, m_narrow, fwhm_narrow)
ax_3d_w = fig.add_subplot(gs[1, 1], projection="3d")
surf2 = style_axis(ax_3d_w, image_wide, m_wide, fwhm_wide)

# Disk/doughnut legend at the BOTTOM of the first row (in the gap between rows)
legend_y = (ax_2d_n.get_position().y0 + ax_3d_n.get_position().y1) / 2
fig.legend(
    handles=[Patch(facecolor="none", edgecolor=DISK_COLOR, lw=2, label="Disk (signal)"),
             Patch(facecolor="none", edgecolor=DONUT_COLOR, lw=2, label="Doughnut (background)")],
    loc="center", bbox_to_anchor=(0.5, legend_y), ncol=2, fontsize=ANNOT_SIZE, frameon=True)

# Colorbar at the very bottom
cbar_ax = fig.add_axes([0.35, 0.055, 0.3, 0.018])
cbar = fig.colorbar(surf1, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Intensity", fontsize=LABEL_SIZE, color="black")
cbar.ax.tick_params(labelsize=TICK_SIZE, colors="black")

output_path_comparison = results_folder / "synthetic_spot_size_vs_intensity.png"
fig.savefig(output_path_comparison, dpi=300, bbox_inches="tight")
plt.show()

# ---------------------------------------------------------------------
# Print numerical values
# ---------------------------------------------------------------------
print("Synthetic translation spot analysis comparison")
print("----------------------------------------------")
print(f"Field of view: {CROP_SIZE} x {CROP_SIZE} px")
print(f"Disk side (spot_size + 2): {best_size} px | doughnut crop: {best_size + PIXELS_AROUND_SPOT} px")
print(f"SNR method: '{SNR_METHOD}' (mean-based); target SNR = {SNR_TARGET}")
print(f"Construction: matched disk-doughnut intensity AND matched SNR; background free\n")

for label, sigma, amp, fwhm, m, nsd in [
    ("Narrow Spot (typical / UTag-like)", SIGMA_NARROW, AMPLITUDE_NARROW, fwhm_narrow, m_narrow, noise_sd_narrow),
    ("Wide Spot (longer array / SunTag-like)", SIGMA_WIDE, AMPLITUDE_WIDE, fwhm_wide, m_wide, noise_sd_wide),
]:
    print(f"{label}:")
    print(f"  Gaussian sigma:               {sigma:.2f} px")
    print(f"  Gaussian peak amplitude:      {amp:.2f}")
    print(f"  Gaussian FWHM (spot size):    {fwhm:.2f} px")
    print(f"  Background-noise sd (free):   {nsd:.2f}")
    print(f"  Disk mean:                    {m['mean_disk']:.2f}")
    print(f"  Doughnut mean (measured bg):  {m['mean_donut']:.2f}")
    print(f"  Doughnut std:                 {m['std_donut']:.2f}")
    print(f"  Spot intensity (disk-donut):  {m['intensity']:.2f}")
    print(f"  SNR (mean / disk-doughnut):   {m['snr_disk_donut']:.2f}")
    print(f"  SNR (peak):                   {m['snr_peak']:.2f}\n")

d_fwhm = (fwhm_wide / fwhm_narrow - 1) * 100
print("Key Observation:")
print(f"  With matched intensity ({m_narrow['intensity']:.1f} vs {m_wide['intensity']:.1f}) and matched SNR "
      f"({m_narrow['snr_disk_donut']:.1f} vs {m_wide['snr_disk_donut']:.1f}),")
print(f"  the spot size (FWHM) is {d_fwhm:+.1f}% larger ({fwhm_narrow:.2f} -> {fwhm_wide:.2f} px) for the")
print(f"  broader spot. Its peak amplitude is LOWER ({AMPLITUDE_WIDE:.0f} vs {AMPLITUDE_NARROW:.0f}) and its")
print(f"  local background differs ({m_narrow['mean_donut']:.1f} vs {m_wide['mean_donut']:.1f}). A larger")
print(f"  apparent spot size therefore does not imply higher fluorescence or SNR.")
