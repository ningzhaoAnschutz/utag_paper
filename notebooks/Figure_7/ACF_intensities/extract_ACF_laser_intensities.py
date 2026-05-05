"""
Extract Laser Intensities from ACF LIF Files

This script processes each subfolder in the ACF data directory
independently, extracting laser intensity metadata from all LIF files and
saving a separate Excel file per folder.

Usage:
    python extract_ACF_laser_intensities.py

Output (one Excel per folder):
    - AlfaTag_laser_intensities.xlsx
    - SunTag_laser_intensities.xlsx
    - UTag_laser_intensities.xlsx
    - UTag_CF_laser_intensities.xlsx
"""

from pathlib import Path
import pandas as pd
from microlive import microscopy as mi


# ── Configuration ──────────────────────────────────────────────────────────────
DATA_ROOT = Path("/Volumes/Luis_DRIVE/UTag_paper_data/ACF")
OUTPUT_DIR = Path(__file__).resolve().parent  # same folder as this script


# ── Core extraction function ──────────────────────────────────────────────────
import sys

# Add repo root to path to access src/utilities
repo_root = Path(__file__).resolve().parent.parent.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.utilities.extract_laser_intensities import extract_laser_intensities


def save_to_excel(df, excel_path):
    """Save DataFrame to Excel with auto-fitted column widths."""
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Laser Intensities", index=False)

        # Auto-fit column widths for readability
        worksheet = writer.sheets["Laser Intensities"]
        for col_idx, col_name in enumerate(df.columns, start=1):
            max_len = max(
                len(str(col_name)),
                df[col_name].astype(str).str.len().max()
                if not df[col_name].isna().all()
                else 0,
            )
            worksheet.column_dimensions[
                worksheet.cell(row=1, column=col_idx).column_letter
            ].width = min(max_len + 3, 50)

    print(f"  ✅ Excel saved to: {excel_path}")


def print_summary(df, folder_name):
    """Print a per-file intensity summary."""
    if df.empty:
        return

    intensity_cols = [c for c in df.columns if "Intensity (%)" in c]

    print(f"\n  {'─' * 60}")
    print(f"  Summary for {folder_name}:")

    for col in intensity_cols:
        summary = df.groupby("File Name").agg(
            images=(col, "count"),
            mean=(col, "mean"),
            min_val=(col, "min"),
            max_val=(col, "max"),
        ).reset_index()

        print(f"    Laser: {col}")
        for _, row in summary.iterrows():
            if row["min_val"] == row["max_val"]:
                print(f"      {row['File Name']}: {row['images']} images, {row['min_val']}%")
            else:
                print(f"      {row['File Name']}: {row['images']} images, {row['min_val']}–{row['max_val']}% (mean {row['mean']:.1f}%)")
    print()


# ── Main execution ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("Laser Intensity Extraction — ACF Dataset")
    print("=" * 70)
    print(f"\nData root  : {DATA_ROOT}")
    print(f"Output dir : {OUTPUT_DIR}\n")

    if not DATA_ROOT.exists():
        raise FileNotFoundError(f"Data root not found:\n  {DATA_ROOT}")

    # Discover subfolders (skip hidden files/dirs)
    subfolders = sorted(
        d for d in DATA_ROOT.iterdir()
        if d.is_dir() and not d.name.startswith(".")
    )

    print(f"Found {len(subfolders)} subfolders: {[d.name for d in subfolders]}\n")

    # Process each folder independently
    for folder in subfolders:
        print("=" * 70)
        print(f"Processing folder: {folder.name}")
        print("=" * 70)

        # Extract intensities
        df = extract_laser_intensities(folder, verbose=False)

        if df.empty:
            continue

        # Save Excel
        excel_path = OUTPUT_DIR / f"{folder.name}_laser_intensities.xlsx"
        save_to_excel(df, excel_path)

        # Print summary
        print_summary(df, folder.name)

    print("\n✅ All folders processed.")
