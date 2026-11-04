"""
Thin DAS read wrapper that also returns the TRUE start time of the returned data.

The vendored `Read_DAS_data` returns only the array; the underlying `load_das_custom`
does not front-pad a request that begins before the earliest available sample (or across a
gap), so the caller cannot assume the array starts exactly at the requested t0. The picker
needs the real start time to place its display window and per-channel arrivals correctly at
file/gap boundaries. This wrapper mirrors `Read_DAS_data`'s slicing + preprocessing exactly
but also hands back `axis["date_times"][0]`.
"""

from readers import preprocess                       # vendored
import load_das_h5_PE2024 as _loader                 # vendored


def read_window(t0, t1, cfg):
    """Return (data [nt, nch] preprocessed, t_actual0 datetime = true start of `data`)."""
    data, _headers, axis = _loader.load_das_custom(t0, t1, input_dir=cfg.das_data_dir)
    data = data[:, int(cfg.d1 / cfg.dx):int(cfg.d2 / cfg.dx)]
    data = preprocess(data, cfg.fs, cfg.f1, cfg.f2)
    t_actual0 = axis["date_times"][0]
    return data, t_actual0
