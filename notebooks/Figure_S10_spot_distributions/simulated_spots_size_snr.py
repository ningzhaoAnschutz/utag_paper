"""Synthetic spot illustration measured with the current microLive pipeline.

The script follows three steps:

1. Generate two Gaussian spots with additive noise.
2. Measure both spots with the current microLive ``Intensity`` class.
3. Validate, plot, and report only measurements returned by
   ``microlive.microscopy.Intensity.calculate_intensity``.

This is an illustration of the measurement definitions, not a mechanistic
model of the SunTag reporter.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from microlive import __version__ as MICROLIVE_VERSION
from microlive import microscopy as mi


# Synthetic-image settings
CROP_SIZE = 13
BASELINE = 100.0
SEED = 64

# Fixed generating parameters selected with microLive 1.0.44. Both spots use
# independent draws from the same noise distribution. The wide spot has a
# lower peak amplitude and a larger width; microLive returns nearly identical
# intensity and disk-doughnut SNR, but an approximately 15% larger FWHM.
NARROW_PARAMETERS = (1.3000, 300.00, 12.00)  # sigma, amplitude, noise SD
WIDE_PARAMETERS = (1.4134, 265.17, 12.00)    # sigma, amplitude, noise SD

# microLive settings used for the measurements
SPOT_SIZE = 5
SNR_METHOD = "disk_doughnut"
FAST_GAUSSIAN_FIT = False

# Output settings
OUTPUT_DIR = Path(__file__).parent / "results_spots_synthetic"
OUTPUT_PATH = OUTPUT_DIR / "synthetic_spot_size_vs_intensity.png"

FWHM_FACTOR = 2 * np.sqrt(2 * np.log(2))
CENTER_INDEX = CROP_SIZE // 2
COORDINATES_ZYX = np.array([[0, CENTER_INDEX, CENTER_INDEX]], dtype=int)

pixel_axis = np.arange(CROP_SIZE)
X, Y = np.meshgrid(pixel_axis, pixel_axis)
R2 = (X - CENTER_INDEX) ** 2 + (Y - CENTER_INDEX) ** 2


def generate_spot(sigma, amplitude, noise_sd, noise_pattern):
    """Return one floating-point Gaussian spot with additive noise."""
    signal = amplitude * np.exp(-R2 / (2 * sigma**2))
    return BASELINE + signal + noise_sd * noise_pattern


def measure_with_microlive(image):
    """Measure one synthetic spot with microLive's current Intensity class."""
    image_zyxc = image[None, :, :, None]
    results = mi.Intensity(
        original_image=image_zyxc,
        spot_size=SPOT_SIZE,
        array_spot_location_z_y_x=COORDINATES_ZYX,
        use_max_projection=False,
        optimize_spot_size=False,
        allow_subpixel_repositioning=False,
        fast_gaussian_fit=FAST_GAUSSIAN_FIT,
        snr_method=SNR_METHOD,
    ).calculate_intensity()

    (
        intensities,
        _,
        snrs,
        background_means,
        background_stds,
        fitted_amplitudes,
        fitted_sigmas,
        _,
    ) = results

    fitted_sigma = float(fitted_sigmas[0, 0])
    return {
        "intensity": float(intensities[0, 0]),
        "snr": float(snrs[0, 0]),
        "background_mean": float(background_means[0, 0]),
        "background_std": float(background_stds[0, 0]),
        "fitted_amplitude": float(fitted_amplitudes[0, 0]),
        "fitted_sigma": fitted_sigma,
        "fwhm": FWHM_FACTOR * fitted_sigma,
    }


def validate_measurements(narrow, wide):
    """Fail clearly if a future microLive change alters the illustration."""
    relative_intensity_difference = abs(
        wide["intensity"] / narrow["intensity"] - 1
    )
    relative_snr_difference = abs(wide["snr"] / narrow["snr"] - 1)
    fwhm_ratio = wide["fwhm"] / narrow["fwhm"]

    if relative_intensity_difference > 0.01:
        raise RuntimeError("Synthetic spot intensities differ by more than 1%.")
    if relative_snr_difference > 0.01:
        raise RuntimeError("Synthetic spot SNRs differ by more than 1%.")
    if not 1.14 <= fwhm_ratio <= 1.16:
        raise RuntimeError("Synthetic spot FWHM difference is no longer about 15%.")


def configure_plot_style():
    """Apply the figure's publication text style."""
    plt.rcParams.update(
        {
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
        }
    )


def plot_top_view(ax, image, label, measurement, global_max):
    """Plot a pixel-resolved top view of one measured synthetic spot."""
    extent = [-0.5, CROP_SIZE - 0.5, -0.5, CROP_SIZE - 0.5]
    shown = ax.imshow(
        image,
        cmap="coolwarm",
        origin="lower",
        extent=extent,
        vmin=0,
        vmax=global_max,
    )
    ax.set_xticks(np.arange(CROP_SIZE))
    ax.set_yticks(np.arange(CROP_SIZE))
    ax.set_xticks(np.arange(-0.5, CROP_SIZE, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, CROP_SIZE, 1), minor=True)
    ax.grid(which="minor", color="0.55", linewidth=0.6)
    ax.tick_params(labelsize=8)
    ax.set_xlim(-0.5, CROP_SIZE - 0.5)
    ax.set_ylim(-0.5, CROP_SIZE - 0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("X (px)", fontsize=10)
    ax.set_ylabel("Y (px)", fontsize=10)
    ax.set_title(
        f"{label} (microLive FWHM = {measurement['fwhm']:.2f} px)",
        fontsize=12,
    )
    return shown


def plot_surface(ax, image, measurement, global_max):
    """Plot one synthetic spot and annotate microLive measurements."""
    surface = ax.plot_surface(
        X,
        Y,
        image,
        cmap="coolwarm",
        edgecolor="none",
        vmin=0,
        vmax=global_max,
    )
    annotation = (
        f"Intensity: {measurement['intensity']:.1f}\n"
        f"SNR: {measurement['snr']:.1f}\n"
        f"FWHM: {measurement['fwhm']:.2f} px\n"
        f"Background: {measurement['background_mean']:.1f}"
    )
    ax.text2D(0.04, 0.80, annotation, transform=ax.transAxes, fontsize=10)
    ax.set_xlim(0, CROP_SIZE - 1)
    ax.set_ylim(0, CROP_SIZE - 1)
    ax.set_zlim(0, global_max)
    ax.set_box_aspect((1, 1, 0.8))
    ax.set_xlabel("X (px)", fontsize=10)
    ax.set_ylabel("Y (px)", fontsize=10)
    ax.set_zlabel("Intensity", fontsize=10)
    ax.tick_params(axis="both", labelsize=8)
    ax.tick_params(axis="z", labelsize=8)
    ax.view_init(elev=25, azim=-45)
    return surface


def save_figure(image_narrow, image_wide, measurement_narrow, measurement_wide):
    """Create and save the two-spot comparison figure."""
    global_max = max(float(image_narrow.max()), float(image_wide.max())) + 50

    figure = plt.figure(figsize=(10, 10))
    grid = GridSpec(
        2,
        2,
        figure=figure,
        height_ratios=[0.85, 1.1],
        hspace=0.20,
        wspace=0.22,
        left=0.08,
        right=0.95,
        top=0.93,
        bottom=0.12,
    )

    top_narrow = figure.add_subplot(grid[0, 0])
    top_wide = figure.add_subplot(grid[0, 1])
    plot_top_view(
        top_narrow,
        image_narrow,
        "Narrow",
        measurement_narrow,
        global_max,
    )
    plot_top_view(
        top_wide,
        image_wide,
        "Wide",
        measurement_wide,
        global_max,
    )

    surface_narrow = figure.add_subplot(grid[1, 0], projection="3d")
    surface_wide = figure.add_subplot(grid[1, 1], projection="3d")
    color_source = plot_surface(
        surface_narrow,
        image_narrow,
        measurement_narrow,
        global_max,
    )
    plot_surface(surface_wide, image_wide, measurement_wide, global_max)

    colorbar_axis = figure.add_axes([0.35, 0.055, 0.3, 0.018])
    colorbar = figure.colorbar(color_source, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("Intensity", fontsize=10)
    colorbar.ax.tick_params(labelsize=8)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight")
    plt.close(figure)


def print_results(
    narrow_parameters,
    wide_parameters,
    narrow_measurement,
    wide_measurement,
):
    """Print the generating parameters and microLive measurements."""
    size_difference = (
        wide_measurement["fwhm"] / narrow_measurement["fwhm"] - 1
    ) * 100

    print("Synthetic translation spots measured with current microLive")
    print("-----------------------------------------------------------")
    print(f"microLive version: {MICROLIVE_VERSION}")
    print(f"microLive SNR method: {SNR_METHOD}")
    print(f"microLive fast_gaussian_fit: {FAST_GAUSSIAN_FIT}")
    print(f"microLive spot_size input: {SPOT_SIZE} px\n")

    for label, parameters, measurement in (
        ("Narrow", narrow_parameters, narrow_measurement),
        ("Wide", wide_parameters, wide_measurement),
    ):
        sigma, amplitude, noise_sd = parameters
        print(f"{label} spot:")
        print(f"  Input sigma:                 {sigma:.4f} px")
        print(f"  Input peak amplitude:        {amplitude:.2f}")
        print(f"  Input noise SD:              {noise_sd:.2f}")
        print(f"  microLive intensity:         {measurement['intensity']:.4f}")
        print(f"  microLive SNR:               {measurement['snr']:.4f}")
        print(f"  microLive background mean:   {measurement['background_mean']:.4f}")
        print(f"  microLive background SD:     {measurement['background_std']:.4f}")
        print(f"  microLive fitted amplitude:  {measurement['fitted_amplitude']:.4f}")
        print(f"  microLive fitted sigma:      {measurement['fitted_sigma']:.4f} px")
        print(f"  microLive FWHM:              {measurement['fwhm']:.4f} px\n")

    print("Key observation:")
    print(
        "  microLive measured nearly identical intensity "
        f"({narrow_measurement['intensity']:.2f} vs "
        f"{wide_measurement['intensity']:.2f}) and SNR "
        f"({narrow_measurement['snr']:.2f} vs {wide_measurement['snr']:.2f})."
    )
    print(
        f"  microLive measured a {size_difference:.1f}% larger FWHM "
        f"({narrow_measurement['fwhm']:.2f} vs "
        f"{wide_measurement['fwhm']:.2f} px) for the wider spot."
    )
    print(f"  Figure saved to: {OUTPUT_PATH}")


def main():
    """Generate, measure, plot, and report the synthetic spots."""
    configure_plot_style()
    rng = np.random.default_rng(SEED)
    noise_narrow = rng.standard_normal((CROP_SIZE, CROP_SIZE))
    noise_wide = rng.standard_normal((CROP_SIZE, CROP_SIZE))

    image_narrow = generate_spot(*NARROW_PARAMETERS, noise_narrow)
    measurement_narrow = measure_with_microlive(image_narrow)

    image_wide = generate_spot(*WIDE_PARAMETERS, noise_wide)
    measurement_wide = measure_with_microlive(image_wide)
    validate_measurements(measurement_narrow, measurement_wide)

    save_figure(
        image_narrow,
        image_wide,
        measurement_narrow,
        measurement_wide,
    )
    print_results(
        NARROW_PARAMETERS,
        WIDE_PARAMETERS,
        measurement_narrow,
        measurement_wide,
    )


if __name__ == "__main__":
    main()
