"""
v3 detection = the v2 detection loops with ONE injection point: after the STA/LTA
characteristic function is computed (and trimmed), the masked channels' CF columns are
set to 0 BEFORE triggering/tracking.

Why this seam (measured on the demo data, see V3_CHANGES.md):
  * zeroing the RAW data column instead -> the MAD-normalised STA/LTA turns the whole
    column into NaN, which pollutes downstream percentiles/plots;
  * a constant/degenerate channel left alone -> inf CF and FALSE triggers;
  * zeroing the CF column -> trig_onechan finds nothing, the channel counts as "empty",
    and method B's chain tracking jumps over it exactly like a dead channel
    (neighbor_tol, default 10). No NaN anywhere, no kernel modification.

Everything else is a faithful copy of dasapp/detect.py (which is NOT modified). The copy
is deliberate: the v2 loop is monolithic (no injection hook), and v3 must not edit it.
Differences from the v2 file are marked with `# V3:`.

The catalog gains two provenance columns (harmless to v2 readers, which select columns
by name): mask_n (number of masked channels) and mask_hash. The mask actually used is
also copied into the experiment folder as channel_mask_used.csv.
"""

import os
import gc
import shutil
import numpy as np
import pandas as pd
from datetime import timedelta

from stalta import (                                    # vendored (v2 path)
    sta_lta_centred_multichan,
    trig_allchan_and_collapsed_index,
    trig_allchan_and_collapsed_index_and_save_wavefonts,
)
from dasapp.detect import _chunks, _flush               # v2 helpers, reused untouched
from .qcio import read_das_qc, apply_csl_mask           # V3: sanitising reader + seam

try:
    import h5py
except Exception:
    h5py = None


def _mask_meta(cfg):
    mask = cfg.chan_mask
    n = 0 if mask is None else int(mask.sum())
    return n, (cfg.mask_tag if n else "nomask")


def _effective_min_channels(cfg, nch):
    """V3: scale min_channels by the surviving-channel fraction (once, globally),
    so the gate keeps the same MEANING when part of the array is masked."""
    if cfg.chan_mask is None or not cfg.quality.get("scale_min_channels", True):
        return int(cfg.min_channels)
    n_bad = int(cfg.chan_mask[:nch].sum())
    return max(2, round(cfg.min_channels * (nch - n_bad) / max(nch, 1)))


def run_detection_v3(cfg, progress=True, resume=True):
    """Masked detection per cfg.method. cfg is a ConfigV3 with chan_mask resolved."""
    os.makedirs(cfg.exp_dir, exist_ok=True)
    cfg.save_snapshot()
    if cfg.chan_mask_path:                              # V3: provenance copy
        try:
            shutil.copy(cfg.chan_mask_path,
                        os.path.join(cfg.exp_dir, "channel_mask_used.csv"))
        except Exception:
            pass
    if cfg.method == "B" and h5py is None:
        raise RuntimeError("h5py is required for method B (per-event wavefront H5).")
    dispatch = {"A": _detect_simple_v3, "B": _detect_wavefront_v3}
    if cfg.method not in dispatch:
        raise ValueError(f"Unknown method {cfg.method!r} (use 'A' or 'B').")
    return dispatch[cfg.method](cfg, progress=progress, resume=resume)


# ---------------------------------------------------------------- Method B ------------
def _detect_wavefront_v3(cfg, progress=True, resume=True):
    os.makedirs(cfg.wavefront_dir, exist_ok=True)
    chunks = _chunks(cfg)
    mask = cfg.chan_mask                                # V3
    mask_n, mask_hash = _mask_meta(cfg)                 # V3

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
        data, dead_raw = read_das_qc(t0_pad, t1, cfg)    # V3: sanitising reader
        csl_pad = sta_lta_centred_multichan(data, cfg.nsta, cfg.nlta)
        csl = csl_pad[cfg.nsta:, :]                      # trim the pad -> starts at t0
        # V3: the one injection point — static mask ∪ this chunk's degenerate channels
        chunk_mask = dead_raw if mask is None else (mask[: len(dead_raw)] | dead_raw)
        csl = apply_csl_mask(csl, chunk_mask)
        mc_eff = _effective_min_channels(cfg, csl.shape[1])          # V3

        (idx_starts, idx_ends, scores, avg_vals, empty_chan,
         wavefronts) = trig_allchan_and_collapsed_index_and_save_wavefonts(
            csl, cfg.threshold_level, cfg.minimal_event_interval, cfg.fs,
            cfg.dx, cfg.tracking_params, cfg.merge_params)

        n_valid = 0
        for wf in wavefronts:
            if wf["score"] < mc_eff:                     # V3: scaled gate
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
                    g.attrs["mask_n"] = mask_n                       # V3
                    g.attrs["mask_hash"] = mask_hash                 # V3
                    g.create_dataset("wavefront_data",
                                     data=np.array(wf["details"]), compression="gzip")
            buffer.append({
                "start_time": abs_start, "end_time": abs_end,
                "duration_s": (abs_end - abs_start).total_seconds(),
                "n_channels": wf["score"], "avg_amp": wf["peak_val_sum"] / wf["score"],
                "method": "B", "h5_file": h5_file, "h5_group": grp,
                "mask_n": mask_n, "mask_hash": mask_hash,            # V3
            })
        if progress:
            n_msk = int(chunk_mask.sum())
            print(f"    -> {n_valid} events (raw {len(wavefronts)}; minCh_eff {mc_eff}; "
                  f"masked {n_msk} ch; empty_unmasked {empty_chan - n_msk})")   # V3

        if buffer and (i % 20 == 0 or i == len(chunks) - 1):
            existing = _flush(cfg.catalog_path, buffer, existing)
            buffer = []
            gc.collect()

    existing = _flush(cfg.catalog_path, buffer, existing)
    n = 0 if isinstance(existing, pd.DataFrame) and existing.empty else len(existing)
    if progress:
        print(f">>> v3 method B done: {n} events -> {cfg.catalog_path}")
    return cfg.catalog_path


# ---------------------------------------------------------------- Method A ------------
def _detect_simple_v3(cfg, progress=True, resume=True):
    chunks = _chunks(cfg)
    mask = cfg.chan_mask                                # V3
    mask_n, mask_hash = _mask_meta(cfg)                 # V3

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
        data, dead_raw = read_das_qc(t0_pad, t1, cfg)    # V3: sanitising reader
        csl_pad = sta_lta_centred_multichan(data, cfg.nsta, cfg.nlta)
        csl = csl_pad[-data.shape[0]:, :]                # faithful to release method A
        # V3: the one injection point — static mask ∪ this chunk's degenerate channels
        chunk_mask = dead_raw if mask is None else (mask[: len(dead_raw)] | dead_raw)
        csl = apply_csl_mask(csl, chunk_mask)
        mc_eff = _effective_min_channels(cfg, csl.shape[1])          # V3

        idx_max, score, csl_av, empty_no = trig_allchan_and_collapsed_index(
            csl, cfg.threshold_level, cfg.minimal_event_interval, cfg.fs)
        if len(idx_max) == 0:
            continue

        keep = np.where(score > mc_eff)[0]               # V3: scaled gate
        idx_max, score, csl_av = idx_max[keep], score[keep], csl_av[keep]
        if progress:
            n_msk = int(chunk_mask.sum())
            print(f"    -> {len(idx_max)} events (>{mc_eff} ch); "
                  f"masked {n_msk} ch; empty_unmasked {empty_no - n_msk}")   # V3

        for j in range(len(idx_max)):
            peak_time = t0 + timedelta(seconds=idx_max[j] / cfg.fs)
            t0s = peak_time - timedelta(seconds=0.5)
            t1s = peak_time + timedelta(seconds=0.5)
            i0 = int((t0s - t0_pad).total_seconds() * cfg.fs)
            i1 = int((t1s - t0_pad).total_seconds() * cfg.fs)
            i0, i1 = max(0, i0), min(data.shape[0], i1)
            # V3: exclude masked channels from the max-amp summary (their raw data is
            # untouched by the CF masking, but they must not contribute to event stats)
            if i1 > i0:
                sub = np.abs(data[i0:i1, :])
                if chunk_mask.any():
                    sub = sub[:, ~chunk_mask[: sub.shape[1]]]
                max_amp_mean = float(np.mean(np.max(sub, axis=0))) if sub.size else np.nan
            else:
                max_amp_mean = np.nan
            buffer.append({
                "start_time": peak_time, "end_time": peak_time, "duration_s": 0.0,
                "n_channels": int(score[j]), "avg_amp": float(csl_av[j]),
                "method": "A", "h5_file": "", "h5_group": "",
                "peak_amp": float(csl_av[j]), "peak_score": int(score[j]),
                "max_amp_mean": max_amp_mean,
                "mask_n": mask_n, "mask_hash": mask_hash,            # V3
            })

        if buffer and (i % 20 == 0 or i == len(chunks) - 1):
            existing = _flush(cfg.catalog_path, buffer, existing)
            buffer = []
            gc.collect()

    existing = _flush(cfg.catalog_path, buffer, existing)
    n = 0 if isinstance(existing, pd.DataFrame) and existing.empty else len(existing)
    if progress:
        print(f">>> v3 method A done: {n} events -> {cfg.catalog_path}")
    return cfg.catalog_path
