"""
DAS data reading + preprocessing wrapper for the STA/LTA package.

This module is shared (library-level) by both detection mainlines
(method_simple_count / method_wavefront) and the QC scripts.

Three fixes applied relative to the original multichan_stalta_v1130/functions/readers.py
(documented in the package design doc, section 5.3):
  (1) Removed the dead import `from functions.compute_correlation import preprocess`
      (that module does not exist; `preprocess` is defined locally below).
  (2) Changed the absolute module import
        `import multichan_stalta.multichan_stalta_v1130.functions.load_das_h5_PE2024 ...`
      to a relative one (`import load_das_h5_PE2024`); the `functions/` directory is
      placed on sys.path by each top-level script (and by the bootstrap below), so the
      bare module name resolves.
  (3) The hard-coded `/mnt/.../DATA_DEC_80Hz_16m_SORTED/` path is hoisted to an explicit,
      visible variable `DAS_DATA_DIR` at the top of this file (NOT an environment variable,
      so the data flow stays visible). Override it directly here, or pass `input_dir=...`
      to `Read_DAS_data`.
"""

# ! ============================================================================
# !  DATA SOURCE PATH (edit me) -- where the raw DAS .h5 files live.
# !  This was previously hard-coded inside Read_DAS_data; it is now an explicit,
# !  visible top-of-file variable. No environment variables are used on purpose.
# !
# !  By default it points at the lightweight demo sample shipped with this package
# !  (IN_DAS_data/, 2024-04-25 15:00-15:30). For a full run on the cluster, set it
# !  to the original archive, e.g.:
# !      DAS_DATA_DIR = '/mnt/E2_SWP-DAS/202404_GREECE/DATA_DEC_80Hz_16m_SORTED/'
# ! ============================================================================
import os as _os

# functions/ -> package root
_PKG_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))

# Default: the bundled demo sample. load_das_custom() globs recursively for files
# matching the requested date, so pointing at IN_DAS_data/ (which contains 20240425/)
# is sufficient for the demo.
DAS_DATA_DIR = _os.path.join(_PKG_ROOT, 'IN_DAS_data') + _os.sep

# For reference (the original archive locations on the acquisition server):
das_h5_mnt_raw = '/mnt/E2_SWP-DAS/202404_GREECE/DATA_RAW_400Hz_8m/'
das_h5_mnt_80Hz = '/mnt/E2_SWP-DAS/202404_GREECE/DATA_DEC_80Hz_16m_SORTED/'


import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import scipy.signal as signal
from scipy.ndimage import uniform_filter1d
import pandas as pd

import sys
import gc
import pickle
import re
import glob
import h5py
import json

# -- Make the bundled functions/ directory importable so the relative imports below
#    resolve even if a caller forgot to set up sys.path.
_FUNC_DIR = _os.path.dirname(_os.path.abspath(__file__))
if _FUNC_DIR not in sys.path:
    sys.path.insert(0, _FUNC_DIR)

from pydas_readers.readers import load_das_h5, write_das_h5
from pydas_readers.util import block_filters, block_cleaning

# !================= Import the custom load function (FIX 2: relative import) =========
import load_das_h5_PE2024 as load_das_h5_PE2024


def load_pkl(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)


def write_pkl(out_obj, output_file):
    with open(output_file, 'wb') as f:  # 'wb' = write binary mode
        pickle.dump(out_obj, f)


def Read_DAS_data(t_start, t_end_extended, fs, f1, f2, dx, d1, d2, input_dir=None):
    """
    Load and preprocess DAS data for a time window and spatial range.

    Parameters
    ----------
    t_start, t_end_extended : datetime
        Time window to load (t_start should already include any STA/LTA padding).
    fs : float          Sampling rate (Hz).
    f1, f2 : float      Bandpass corners (Hz).
    dx : float          Channel spacing (m).
    d1, d2 : float      Spatial slice start / end distance (m).
    input_dir : str or None
        Directory of raw DAS .h5 files. Defaults to the module-level DAS_DATA_DIR.

    Returns
    -------
    data : 2D np.ndarray, shape (n_time, n_channels), preprocessed.
    """
    if input_dir is None:
        input_dir = DAS_DATA_DIR

    # Load raw data via the custom loader
    data, headers, axis = load_das_h5_PE2024.load_das_custom(
        t_start, t_end_extended, input_dir=input_dir
    )

    # Spatial slicing
    data = data[:, int(d1 / dx):int(d2 / dx)]

    # Preprocessing (spatial median removal + detrend + bandpass)
    data = preprocess(data, fs, f1, f2)
    return data


def preprocess(data, fs, f1, f2):
    """
    INPUT: 2D numpy array, data(i,j): i is time, j is channel number.
    OUTPUT: 2D numpy array, data_pre(i,j), pre-processed, same shape.

    Steps:
    1) Remove the spatial median (common-mode noise reduction across channels).
    2) Remove the linear trend and temporal mean.
    3) Apply a bandpass filter.
    """
    # Remove the spatial median (removes common noise across channels)
    data = data - np.median(data, axis=1)[:, np.newaxis]

    # Remove the linear trend and temporal mean
    data = block_cleaning.detrend(data, type='linear')
    data = data - np.mean(data, axis=0)

    # Apply bandpass filter
    data = block_filters.block_bandpass(data, f1, f2, fs, zerophase=False, taper=0.02)
    return data
