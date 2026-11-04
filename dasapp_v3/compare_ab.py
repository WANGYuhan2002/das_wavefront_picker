"""
A/B comparison of two detection catalogs (typically v2 vs v3 on the same window).

    python run_v3.py compare --a <v2_exp>/catalog.csv --b <v3_exp>/catalog.csv

Reports: kept / lost / gained events, per-event channel-count deltas, and (for method-B
catalogs with Wavefronts/*.h5) tracked-wavefront roughness — the honest way to decide
whether a mask helps on YOUR data before trusting it. Writes ab_report.csv next to
catalog B.
"""

import os
import numpy as np
import pandas as pd

try:
    import h5py
except Exception:
    h5py = None


def _roughness(exp_dir, row, fs):
    """Median |arrival - 7pt median| (samples) of one tracked wavefront, or NaN."""
    if h5py is None or not row.get("h5_file"):
        return np.nan
    from scipy.ndimage import median_filter
    path = os.path.join(exp_dir, "Wavefronts", str(row["h5_file"]))
    try:
        with h5py.File(path, "r") as f:
            det = np.array(f[str(row["h5_group"])]["wavefront_data"])
    except Exception:
        return np.nan
    tt = det[np.argsort(det[:, 0]), 1].astype(float)
    if len(tt) < 9:
        return np.nan
    return float(np.median(np.abs(tt - median_filter(tt, 7, mode="nearest"))))


def _load(cat_path):
    """Load any catalog v2's reader understands (unified, release v1/v2 schemas)."""
    from dasapp import catalog as cat                 # v2, untouched
    events, _schema = cat.load_catalog(cat_path)
    rows = [{"start_time": e["start_time"], "n_channels": e.get("n_channels", 0),
             "h5_file": e.get("h5_file", ""), "h5_group": e.get("h5_group", "")}
            for e in events]
    return pd.DataFrame(rows, columns=["start_time", "n_channels", "h5_file", "h5_group"])


def compare_catalogs(cat_a, cat_b, tol_s=1.0, fs=80.0, label_a="A", label_b="B"):
    """Match events across two catalogs by start_time; return (summary dict, report df)."""
    A = _load(cat_a)
    B = _load(cat_b)
    dir_a, dir_b = (os.path.dirname(os.path.abspath(p)) for p in (cat_a, cat_b))

    used = set()
    rows = []
    for i, a in A.iterrows():
        best_dt, bj = None, None
        for j, b in B.iterrows():
            if j in used:
                continue
            dt = abs((a["start_time"] - b["start_time"]).total_seconds())
            if dt <= tol_s and (best_dt is None or dt < best_dt):
                best_dt, bj = dt, j
        if bj is not None:
            used.add(bj)
            b = B.loc[bj]
            rows.append({
                "start_time": a["start_time"], "status": "matched",
                f"n_ch_{label_a}": a["n_channels"], f"n_ch_{label_b}": b["n_channels"],
                "d_n_ch": int(b["n_channels"] - a["n_channels"]),
                f"rough_{label_a}": _roughness(dir_a, a, fs),
                f"rough_{label_b}": _roughness(dir_b, b, fs),
            })
        else:
            rows.append({"start_time": a["start_time"], "status": f"only_{label_a}",
                         f"n_ch_{label_a}": a["n_channels"],
                         f"rough_{label_a}": _roughness(dir_a, a, fs)})
    for j, b in B.iterrows():
        if j not in used:
            rows.append({"start_time": b["start_time"], "status": f"only_{label_b}",
                         f"n_ch_{label_b}": b["n_channels"],
                         f"rough_{label_b}": _roughness(dir_b, b, fs)})
    if rows:
        rep = pd.DataFrame(rows).sort_values("start_time").reset_index(drop=True)
    else:                                             # both catalogs empty
        rep = pd.DataFrame(columns=["start_time", "status", "d_n_ch",
                                    f"n_ch_{label_a}", f"n_ch_{label_b}",
                                    f"rough_{label_a}", f"rough_{label_b}"])

    matched = rep[rep.status == "matched"]
    ra = matched.get(f"rough_{label_a}", pd.Series(dtype=float)).astype(float)
    rb = matched.get(f"rough_{label_b}", pd.Series(dtype=float)).astype(float)
    summary = {
        f"events_{label_a}": len(A), f"events_{label_b}": len(B),
        "matched": int(len(matched)),
        f"lost (only {label_a})": int((rep.status == f"only_{label_a}").sum()),
        f"gained (only {label_b})": int((rep.status == f"only_{label_b}").sum()),
        "median d_n_ch": float(matched["d_n_ch"].median()) if len(matched) else np.nan,
        f"median rough_{label_a} (samples)": float(np.nanmedian(ra)) if len(ra) else np.nan,
        f"median rough_{label_b} (samples)": float(np.nanmedian(rb)) if len(rb) else np.nan,
        f"p90 rough_{label_a}": float(np.nanpercentile(ra, 90)) if len(ra) else np.nan,
        f"p90 rough_{label_b}": float(np.nanpercentile(rb, 90)) if len(rb) else np.nan,
    }
    return summary, rep


def run_compare(cat_a, cat_b, tol_s=1.0, fs=80.0, out_csv=None):
    summary, rep = compare_catalogs(cat_a, cat_b, tol_s=tol_s, fs=fs,
                                    label_a="v2", label_b="v3")
    out_csv = out_csv or os.path.join(os.path.dirname(os.path.abspath(cat_b)),
                                      "ab_report.csv")
    rep.to_csv(out_csv, index=False)
    print("=== A/B comparison ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    lost = rep[rep.status == "only_v2"]
    if len(lost):
        print("  LOST events (investigate before trusting the mask):")
        for _, r in lost.iterrows():
            print(f"    {r['start_time']}  n_ch={r.get('n_ch_v2')}")
    gained = rep[rep.status == "only_v3"]
    if len(gained):
        print("  gained events:")
        for _, r in gained.iterrows():
            print(f"    {r['start_time']}  n_ch={r.get('n_ch_v3')}")
    print(f"  report -> {out_csv}")
    return summary, out_csv
