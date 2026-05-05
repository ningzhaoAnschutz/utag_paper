"""Shared imports, path setup, and colormap definitions for utag_paper notebooks."""

# Standard library
import sys
import os
import pathlib
from pathlib import Path
import re
import math
import time
import random
import datetime
from datetime import datetime
from itertools import combinations
import subprocess
import warnings
import logging

# Scientific computing
import numpy as np
import pandas as pd
import cv2
from scipy.optimize import curve_fit
from scipy.stats import linregress, pearsonr, mannwhitneyu
import scipy.stats as stats
from scipy.ndimage import gaussian_filter
from joblib import Parallel, delayed

# Visualization
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_scalebar.scalebar import ScaleBar
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from PIL import Image

# Microscopy I/O
import tifffile
from skimage import exposure
from skimage.measure import label

# Jupyter
from IPython.display import display
import ipywidgets as widgets

# Cellpose (suppress verbose output)
import contextlib, io
logging.getLogger("root").setLevel(logging.ERROR)
_f = io.StringIO()
with contextlib.redirect_stdout(_f), contextlib.redirect_stderr(_f):
    from cellpose import models

# Microlive — import from the pip-installed package
try:
    import microlive.microscopy as mi
except ImportError:
    try:
        # Fallback: legacy path-based import
        _microlive_dir = Path.home() / 'Desktop' / 'microlive' / 'src'
        if _microlive_dir.is_dir():
            sys.path.append(str(_microlive_dir))
        import microscopy as mi
    except ImportError:
        print("Warning: microlive not available. Some functionality will be limited.")

# --- Project path setup ---
warnings.filterwarnings("ignore")

def _find_project_root():
    """Walk up from cwd to find the utag project root (contains src/ and notebooks/)."""
    check = Path.cwd()
    while check != check.parent:
        if check.name in ('utag', 'utag_paper') and (check / 'src').exists():
            return check
        if (check / 'src').exists() and (check / 'notebooks').exists():
            return check
        check = check.parent
    # Fallback locations
    for loc in [Path.home() / 'Desktop' / 'utag_paper', Path.home() / 'Desktop' / 'utag']:
        if loc.exists() and (loc / 'src').exists():
            return loc
    return Path.cwd()

main_dir = _find_project_root()
src_dir = main_dir / 'src'
if src_dir.is_dir() and str(src_dir) not in sys.path:
    sys.path.append(str(src_dir))

# --- Custom colormaps (ImageJ-style) ---
def _make_cmap(name, r, g, b):
    return LinearSegmentedColormap(name, {
        'red':   ((0, 0, 0), (1, r, r)),
        'green': ((0, 0, 0), (1, g, g)),
        'blue':  ((0, 0, 0), (1, b, b)),
    })

green_colormap   = _make_cmap('BlackGreen',   0, 1, 0)
magenta_colormap = _make_cmap('BlackMagenta', 1, 0, 1)
red_colormap     = _make_cmap('BlackRed',     1, 0, 0)
yellow_colormap  = _make_cmap('BlackYellow',  1, 1, 0)

cmap_list_imagej = [green_colormap, magenta_colormap, yellow_colormap, red_colormap]