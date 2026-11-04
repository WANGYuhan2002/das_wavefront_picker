"""
Catalog + wavefront-H5 compatibility layer.

Reads *any* STA/LTA catalog produced by the matched_filter_release package or by this
app, auto-detecting the schema:

  * v2 (method B / wavefront):  start_time, end_time, duration_s, n_channels, avg_amp,
                                h5_file, h5_group      (+ per-event Wavefronts/*.h5)
  * v1 (method A / simple):     peak_time, peak_amp, peak_score, max_amp_mean   (no H5)
  * app unified:                same as v2 + a `method` column

Every row is normalised into a plain dict `event` with a stable set of keys so the rest
of the app never has to care which method produced it:

    idx, start_time, end_time, peak_time, n_channels, amp, h5_file, h5_group, source

The per-event wavefront H5 keeps the *exact* release format
(`wavefront_data` = [channel_idx, time_idx_samples, value]) so files written here are a
drop-in for the release's `2_refine_wavefronts.py` / `1b_waterfall_QC.py`.
"""

import os
import re
import glob
import json
import numpy as np
import pandas as pd
import h5py


def _safe_int(x, default=0):
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return default
        return int(x)
    except (ValueError, TypeError):
        return default


def _safe_float(x, default=np.nan):
    try:
        return float(x)
    except (ValueError, TypeError):
        return default


# ------------------------------------------------------------------ load a catalog ----
def load_catalog(path):
    """Return (events, schema) where events is a list of normalised dicts.

    schema is 'v2', 'v1', or 'unified'. Raises FileNotFoundError if missing.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    cols = set(df.columns)

    if {"start_time", "n_channels"} <= cols:
        schema = "unified" if "method" in cols else "v2"
        df["start_time"] = pd.to_datetime(df["start_time"])
        end = pd.to_datetime(df["end_time"], errors="coerce") if "end_time" in cols else df["start_time"]
        end = end.fillna(df["start_time"])              # blank/NaN end_time -> event start
        df = df[df["start_time"].notna()].reset_index(drop=True)
        end = end.reset_index(drop=True)
        events = []
        for i, row in df.iterrows():
            events.append({
                "idx": int(i),
                "start_time": row["start_time"].to_pydatetime(),
                "end_time": end.iloc[i].to_pydatetime(),
                "peak_time": row["start_time"].to_pydatetime(),
                "n_channels": _safe_int(row["n_channels"]),
                "amp": _safe_float(row.get("avg_amp", np.nan)),
                "h5_file": str(row["h5_file"]) if "h5_file" in cols and pd.notna(row.get("h5_file")) else "",
                "h5_group": str(row["h5_group"]) if "h5_group" in cols and pd.notna(row.get("h5_group")) else "",
                "source": str(row.get("method", "B")),
            })
        return events, schema

    if "peak_time" in cols:
        df["peak_time"] = pd.to_datetime(df["peak_time"], errors="coerce")
        df = df[df["peak_time"].notna()].reset_index(drop=True)
        events = []
        for i, row in df.iterrows():
            pk = row["peak_time"].to_pydatetime()
            events.append({
                "idx": int(i),
                "start_time": pk,
                "end_time": pk,
                "peak_time": pk,
                "n_channels": _safe_int(row.get("peak_score", 0)),
                "amp": _safe_float(row.get("peak_amp", np.nan)),
                "h5_file": "",
                "h5_group": "",
                "source": "A",
            })
        return events, "v1"

    raise ValueError(f"Unrecognised catalog schema in {path}: columns={sorted(cols)}")


# ------------------------------------------------- read a method-B seed wavefront ------
def read_seed_wavefront(exp_dir, event):
    """Return the auto-tracked wavefront for `event` as an (N,3) array
    [channel_idx, time_idx_samples, value], or None if no H5 is available.

    `exp_dir` is the experiment folder that holds Wavefronts/<h5_file>.
    """
    h5_file = event.get("h5_file", "")
    grp = event.get("h5_group", "")
    if not h5_file or not grp:
        return None
    # search a couple of plausible locations for the H5
    cands = [
        os.path.join(exp_dir, "Wavefronts", h5_file),
        os.path.join(exp_dir, h5_file),
    ]
    for p in cands:
        if os.path.exists(p):
            try:
                with h5py.File(p, "r") as f:
                    if grp in f and "wavefront_data" in f[grp]:
                        return np.asarray(f[grp]["wavefront_data"][:], dtype=float)
            except Exception:
                return None
    return None


def seed_to_relative(wf, fs):
    """Convert a method-B wavefront [ch, t_idx_samples, val] into per-channel arrivals
    RELATIVE to the event start (seconds): rel = (t_idx - min t_idx)/fs.

    Returns (channels, arrivals_rel_s, values) as float arrays sorted by channel.
    """
    if wf is None or len(wf) == 0:
        return None
    ch = wf[:, 0].astype(float)
    t = wf[:, 1].astype(float)
    val = wf[:, 2].astype(float)
    rel = (t - t.min()) / float(fs)
    order = np.argsort(ch)
    return ch[order], rel[order], val[order]


# ---------------------------------------------------------- write refined output ------
def _safe_label(label):
    """Filesystem-safe label token for the filename."""
    tok = re.sub(r"[^A-Za-z0-9]+", "-", str(label)).strip("-")
    return tok or "unknown"


def write_wavefront_h5(picks_dir, ref_time, label, channels, arrivals_rel_s, values, fs,
                       conf=None, meta=None):
    """Persist ONE wavefront to its OWN file `picks_dir/wf_<ts>_<label>.h5`.

    The file is keyed by the candidate reference time (one window = one wavefront), so
    re-saving after a relabel replaces cleanly: any prior `wf_<ts>_*.h5` for the same
    reference time is removed first (the label lives in the filename). Contents:
      - `wavefront_data`  [channel_idx, time_idx_rel_samples, value]   (release-compatible)
      - `arrival_rel_s`   per-channel arrival in seconds (unambiguous float)
      - `channels`        per-channel index
      - `conf`            per-channel confidence bool (optional)
      - attrs: label, ref_time, fs, n_channels, + anything in `meta`

    Returns the h5 filename (basename).
    """
    os.makedirs(picks_dir, exist_ok=True)
    # ref_time (= the candidate's start_time) uniquely identifies the wavefront: the detector
    # de-duplicates catalog start_times, so no two candidates in one folder share this key.
    ts = ref_time.strftime("%Y%m%d_%H%M%S_%f")
    # remove any earlier file/png for this same wavefront (e.g. a previous label). glob.escape
    # the directory so a metacharacter in the path (e.g. "[") can't defeat the clean replace.
    safe_dir = glob.escape(picks_dir)
    for old in glob.glob(os.path.join(safe_dir, f"wf_{ts}_*.h5")) + \
            glob.glob(os.path.join(safe_dir, f"wf_{ts}_*.png")):
        try:
            os.remove(old)
        except OSError:
            pass

    fname = f"wf_{ts}_{_safe_label(label)}.h5"
    path = os.path.join(picks_dir, fname)

    ch = np.asarray(channels, float)
    rel = np.asarray(arrivals_rel_s, float)
    val = np.asarray(values, float)
    t_idx_rel = np.round(rel * fs).astype(int)
    wf_data = np.column_stack([ch, t_idx_rel, val]).astype(float)

    with h5py.File(path, "w") as f:
        f.attrs["label"] = str(label)
        f.attrs["ref_time"] = ref_time.isoformat()
        f.attrs["fs"] = float(fs)
        f.attrs["n_channels"] = int(len(ch))
        if meta:
            for k, v in meta.items():
                try:
                    f.attrs[k] = v
                except Exception:
                    f.attrs[k] = str(v)
        f.create_dataset("wavefront_data", data=wf_data, compression="gzip")
        f.create_dataset("arrival_rel_s", data=rel, compression="gzip")
        f.create_dataset("channels", data=ch, compression="gzip")
        if conf is not None:
            f.create_dataset("conf", data=np.asarray(conf, bool), compression="gzip")
    return fname


def update_picks_index(path, row):
    """Append/replace one row (one WAVEFRONT) in the picks index CSV, keyed by ref_time,
    so downstream filtering is a one-liner (e.g. df[df.label == 'explosive'])."""
    if os.path.exists(path):
        try:
            df = pd.read_csv(path)
        except Exception:
            df = pd.DataFrame()
        if "ref_time" in df.columns:                 # replace the row for this wavefront
            df = df[df["ref_time"] != row["ref_time"]]
        df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    else:
        df = pd.DataFrame([row])
    df = df.sort_values("ref_time").reset_index(drop=True)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, index=False)
    return path
