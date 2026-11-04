"""
STA/LTA detection orchestration (self-contained; uses the vendored kernel + reader).

Two methods, faithful to the matched_filter_release main loops so counts reproduce:

  * Method A -- simple cross-channel count  (trig_allchan_and_collapsed_index)
        -> events are single peaks: no per-channel wavefront to seed the picker, so the
           picker will auto-pick per channel inside each event window.
  * Method B -- wavefront chain tracking    (..._and_save_wavefonts)
        -> events carry a per-channel wavefront (Wavefronts/Events_<day>.h5, release
           format), which seeds the interactive picker directly.  [default & recommended]

Both write ONE unified catalog.csv:
    start_time, end_time, duration_s, n_channels, avg_amp, method, h5_file, h5_group
    (+ method-A rows also carry peak_amp / peak_score / max_amp_mean for MF compatibility)
"""

import os
import gc
import numpy as np
import pandas as pd
from datetime import timedelta

from readers import Read_DAS_data                       # vendored
from stalta import (                                    # vendored
    sta_lta_centred_multichan,
    trig_allchan_and_collapsed_index,
    trig_allchan_and_collapsed_index_and_save_wavefonts,
)

try:
    import h5py
except Exception:
    h5py = None


def _chunks(cfg):
    t_array = np.arange(cfg.t_start, cfg.t_end, cfg.dt_chunk).astype(object)
    t_array = np.append(t_array, cfg.t_end)
    return [(t_array[i], t_array[i + 1]) for i in range(len(t_array) - 1)]


def _flush(path, rows, existing):
    if not rows:
        return existing
    df_new = pd.DataFrame(rows)
    df = pd.concat([existing, df_new], ignore_index=True) if len(existing) else df_new
    df["start_time"] = pd.to_datetime(df["start_time"])
    df = df.drop_duplicates(subset=["start_time"], keep="last").sort_values("start_time")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)
    return df


def run_detection(cfg, progress=True, resume=True):
    """Run detection per cfg.method. Returns the path to catalog.csv."""
    os.makedirs(cfg.exp_dir, exist_ok=True)
    cfg.save_snapshot()
    if cfg.method == "B" and h5py is None:
        raise RuntimeError("h5py is required for method B (per-event wavefront H5).")

    dispatch = {"A": _detect_simple, "B": _detect_wavefront}
    if cfg.method not in dispatch:
        raise ValueError(f"Unknown method {cfg.method!r} (use 'A' or 'B').")
    return dispatch[cfg.method](cfg, progress=progress, resume=resume)


# ---------------------------------------------------------------- Method B ------------
def _detect_wavefront(cfg, progress=True, resume=True):
    os.makedirs(cfg.wavefront_dir, exist_ok=True)
    chunks = _chunks(cfg)

    existing = pd.DataFrame()
    resume_after = None
    if resume and os.path.exists(cfg.catalog_path):
        existing = pd.read_csv(cfg.catalog_path)
        if not existing.empty and "start_time" in existing:
            existing["start_time"] = pd.to_datetime(existing["start_time"])
            resume_after = existing["start_time"].max()
            if progress:
                print(f">>> resuming after {resume_after}")

    buffer = []
    for i, (t0, t1) in enumerate(chunks):
        if resume_after is not None and t1 < resume_after:
            continue
        if progress:
            print(f"[{i+1}/{len(chunks)}] {t0} -> {t1}")

        t0_pad = t0 - timedelta(seconds=cfg.sta_len)
        data = Read_DAS_data(t0_pad, t1, cfg.fs, cfg.f1, cfg.f2, cfg.dx,
                             cfg.d1, cfg.d2, input_dir=cfg.das_data_dir)
        csl_pad = sta_lta_centred_multichan(data, cfg.nsta, cfg.nlta)
        csl = csl_pad[cfg.nsta:, :]                      # trim the pad -> starts at t0

        (idx_starts, idx_ends, scores, avg_vals, empty_chan,
         wavefronts) = trig_allchan_and_collapsed_index_and_save_wavefonts(
            csl, cfg.threshold_level, cfg.minimal_event_interval, cfg.fs,
            cfg.dx, cfg.tracking_params, cfg.merge_params)

        n_valid = 0
        for wf in wavefronts:
            if wf["score"] < cfg.min_channels:
                continue
            n_valid += 1
            abs_start = t0 + timedelta(seconds=wf["start_index"] / cfg.fs)
            abs_end = t0 + timedelta(seconds=wf["end_index"] / cfg.fs)
            day = abs_start.strftime("%Y%m%d")
            h5_file = f"Events_{day}.h5"
            grp = abs_start.strftime("%H%M%S_%f")
            with h5py.File(os.path.join(cfg.wavefront_dir, h5_file), "a") as f:
                if grp not in f:
                    g = f.create_group(grp)
                    g.attrs["start_time"] = abs_start.isoformat()
                    g.attrs["end_time"] = abs_end.isoformat()
                    g.attrs["duration_sec"] = (abs_end - abs_start).total_seconds()
                    g.attrs["score"] = wf["score"]
                    g.attrs["avg_val"] = wf["peak_val_sum"] / wf["score"]
                    g.create_dataset("wavefront_data",
                                     data=np.array(wf["details"]), compression="gzip")
            buffer.append({
                "start_time": abs_start, "end_time": abs_end,
                "duration_s": (abs_end - abs_start).total_seconds(),
                "n_channels": wf["score"], "avg_amp": wf["peak_val_sum"] / wf["score"],
                "method": "B", "h5_file": h5_file, "h5_group": grp,
            })
        if progress:
            print(f"    -> {n_valid} events (raw {len(wavefronts)}, empty ch {empty_chan})")

        if buffer and (i % 20 == 0 or i == len(chunks) - 1):
            existing = _flush(cfg.catalog_path, buffer, existing)
            buffer = []
            gc.collect()

    existing = _flush(cfg.catalog_path, buffer, existing)
    n = 0 if isinstance(existing, pd.DataFrame) and existing.empty else len(existing)
    if progress:
        print(f">>> method B done: {n} events -> {cfg.catalog_path}")
    return cfg.catalog_path


# ---------------------------------------------------------------- Method A ------------
def _detect_simple(cfg, progress=True, resume=True):
    chunks = _chunks(cfg)

    existing = pd.DataFrame()
    resume_after = None
    if resume and os.path.exists(cfg.catalog_path):
        existing = pd.read_csv(cfg.catalog_path)
        if not existing.empty and "start_time" in existing:
            existing["start_time"] = pd.to_datetime(existing["start_time"])
            resume_after = existing["start_time"].max()

    buffer = []
    for i, (t0, t1) in enumerate(chunks):
        if resume_after is not None and t1 < resume_after:
            continue
        if progress:
            print(f"[{i+1}/{len(chunks)}] {t0} -> {t1}")

        t0_pad = t0 - timedelta(seconds=cfg.sta_len)
        data = Read_DAS_data(t0_pad, t1, cfg.fs, cfg.f1, cfg.f2, cfg.dx,
                             cfg.d1, cfg.d2, input_dir=cfg.das_data_dir)
        csl_pad = sta_lta_centred_multichan(data, cfg.nsta, cfg.nlta)
        csl = csl_pad[-data.shape[0]:, :]                # faithful to release method A

        idx_max, score, csl_av, empty_no = trig_allchan_and_collapsed_index(
            csl, cfg.threshold_level, cfg.minimal_event_interval, cfg.fs)
        if len(idx_max) == 0:
            continue

        keep = np.where(score > cfg.min_channels)[0]
        idx_max, score, csl_av = idx_max[keep], score[keep], csl_av[keep]
        if progress:
            print(f"    -> {len(idx_max)} events (>{cfg.min_channels} ch), empty {empty_no}")

        for j in range(len(idx_max)):
            peak_time = t0 + timedelta(seconds=idx_max[j] / cfg.fs)
            # max amplitude in +-0.5 s window (v1 max_amp_mean, for MF compatibility)
            t0s = peak_time - timedelta(seconds=0.5)
            t1s = peak_time + timedelta(seconds=0.5)
            i0 = int((t0s - t0_pad).total_seconds() * cfg.fs)
            i1 = int((t1s - t0_pad).total_seconds() * cfg.fs)
            i0, i1 = max(0, i0), min(data.shape[0], i1)
            max_amp_mean = float(np.mean(np.max(np.abs(data[i0:i1, :]), axis=0))) if i1 > i0 else np.nan
            buffer.append({
                "start_time": peak_time, "end_time": peak_time, "duration_s": 0.0,
                "n_channels": int(score[j]), "avg_amp": float(csl_av[j]),
                "method": "A", "h5_file": "", "h5_group": "",
                "peak_amp": float(csl_av[j]), "peak_score": int(score[j]),
                "max_amp_mean": max_amp_mean,
            })

        if buffer and (i % 20 == 0 or i == len(chunks) - 1):
            existing = _flush(cfg.catalog_path, buffer, existing)
            buffer = []
            gc.collect()

    existing = _flush(cfg.catalog_path, buffer, existing)
    n = 0 if isinstance(existing, pd.DataFrame) and existing.empty else len(existing)
    if progress:
        print(f">>> method A done: {n} events -> {cfg.catalog_path}")
    return cfg.catalog_path
