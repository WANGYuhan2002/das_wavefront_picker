"""
Raw-domain sanitising DAS reader for v3 (scan + detection).

Why it exists: v2's `preprocess` starts by removing the SPATIAL MEDIAN across channels.
One NaN sample in ONE channel therefore makes that time-row's median NaN and poisons the
whole array for that chunk (then detrend/bandpass smear it everywhere). A constant
(stuck) channel likewise survives preprocessing and produces +inf STA/LTA CF with a
false trigger (measured — see V3_CHANGES.md). v2 has no hook before `preprocess`, so
this module carries a pinned copy of the 3-line `Read_DAS_data` body (readers.py:85-117)
with the sanitiser inserted between slicing and preprocessing. `_pins.py` hash-guards
the copied original so future v2 edits surface loudly.

With clean data (no NaN, no stuck channels) the output is bit-identical to v2's
`Read_DAS_data` — the maskless v3 detection path stays exactly reproducible.
"""

import numpy as np

import load_das_h5_PE2024 as _loader                  # vendored (v2 path)
from readers import preprocess                        # vendored (v2 path)


def read_das_qc(t_start, t_end, cfg):
    """Load + sanitise + preprocess one window.

    Returns (data [nt, nch] preprocessed, dead_raw bool[nch]) where dead_raw flags
    channels that were non-finite or constant in the RAW chunk — such channels must be
    added to the chunk's CF mask (their content after zero-filling is meaningless).
    """
    data, _headers, _axis = _loader.load_das_custom(t_start, t_end,
                                                    input_dir=cfg.das_data_dir)
    data = data[:, int(cfg.d1 / cfg.dx):int(cfg.d2 / cfg.dx)]

    nonfin = ~np.isfinite(data).all(axis=0)
    ptp0 = (np.ptp(data, axis=0) == 0)
    dead_raw = nonfin | ptp0
    if nonfin.any():
        data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)

    data = preprocess(data, cfg.fs, cfg.f1, cfg.f2)
    return data, dead_raw


def apply_csl_mask(csl, mask):
    """THE injection seam: zero the CF columns of masked channels (in place).

    Guarantees (measured): csl >= 0 with threshold > 0 -> no triggers -> the channel
    counts as `empty_chan` and method-B tracking bridges it via neighbor_tol, exactly
    like a genuinely silent channel. No NaN/inf is ever introduced.
    """
    if mask is not None and mask.any():
        csl[:, mask[: csl.shape[1]]] = 0.0
    return csl
