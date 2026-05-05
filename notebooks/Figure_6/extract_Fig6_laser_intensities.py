"""
Extract Laser Intensities from LIF Files — Anti-Utag Frankenbody Translocation Dataset

This script extracts laser intensity metadata from all LIF files in the specified
directory and saves the results as an Excel file. It dynamically imports the core
extraction logic from the unified utilities module.

Usage:
    python extract_Fig6_laser_intensities.py

Output:
    - laser_intensities.xlsx  (saved in the same directory as this script)
"""

import sys
from pathlib import Path
import pandas as pd

# Add repo root to path to access src/utilities
repo_root = Path(__file__).resolve().parent.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from src.utilities.extract_laser_intensities import extract_laser_intensities


# ── Configuration ──────────────────────────────────────────────────────────────
DATA_DIR = Path(
    "/Users/nzlab-la/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/"
    "General - Zhao (NZ) Lab/Zhao lab shared folder/Our papers/"
    "Anti-Utag-frankenbody paper/Fig. S Translocation/final"
)

OUTPUT_DIR = Path(__file__).resolve().parent  # same folder as this script


# ── Main execution ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("Laser Intensity Extraction — Anti-Utag Translocation Dataset")
    print("=" * 70)
    print(f"\nData directory : {DATA_DIR}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    if not DATA_DIR.exists():
        raise FileNotFoundError(f"Data directory not found:\n  {DATA_DIR}")

    # Extract intensities
    df_intensities = extract_laser_intensities(DATA_DIR, verbose=True)

    # ── Save to Excel ──────────────────────────────────────────────────────
    excel_path = OUTPUT_DIR / "laser_intensities.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        df_intensities.to_excel(writer, sheet_name="Laser Intensities", index=False)

        # Auto-fit column widths for readability
        worksheet = writer.sheets["Laser Intensities"]
        for col_idx, col_name in enumerate(df_intensities.columns, start=1):
            max_len = max(
                len(str(col_name)),
                df_intensities[col_name].astype(str).str.len().max()
                if not df_intensities[col_name].isna().all()
                else 0,
            )
            worksheet.column_dimensions[
                worksheet.cell(row=1, column=col_idx).column_letter
            ].width = min(max_len + 3, 50)

    print(f"\n✅ Excel saved to: {excel_path}")

    # Also save CSV as a backup
    csv_path = OUTPUT_DIR / "laser_intensities.csv"
    df_intensities.to_csv(csv_path, index=False)
    print(f"✅ CSV saved to:   {csv_path}")

    # Quick summary
    if not df_intensities.empty:
        print(f"\n{'─' * 70}")
        print(f"Summary: {len(df_intensities)} total images from {DATA_DIR.name}")
        print(f"Columns: {list(df_intensities.columns)}")
        print(df_intensities.head(10).to_string(index=False))
