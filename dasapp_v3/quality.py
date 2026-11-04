"""
Channel-quality scan: measure per-channel health from continuous DAS data and derive a
bad-channel mask for detection/picking. This module is the heart of v3.

Two distinct notions of channel quality are handled — they need DIFFERENT statistics.
Every rule below was validated (or overturned) empirically on the bundled Kefalonia
demo; the measurements are in V3_CHANGES.md.

  STATIC quality — is the channel healthy in general?  (BACKGROUND windows only)
      * dead  : non-finite / constant raw content, or background z <= z_dead
      * noisy : background z >= z_noisy, or the corroborated soft rule
                (z >= z_noisy_soft AND trigger-rate z >= z_rate AND rate >= rate_min);
                validated: the persistently-hot 13.5-13.9 km cluster sits at z +3..+5
                in every background window
      * trigger rate alone is REPORT-ONLY: on active windows a high trigger rate means
        a GOOD channel (it participates in real events) — an absolute threshold masked
        30% of the demo array, including its best channels.

  RESPONSIVENESS — does the channel see energy when the array does? (EXCITED windows)
      * weak responder: global-z < resp_z in >= resp_min_votes excited windows.
        Measured facts that fixed this rule's shape:
          - weak channels are INVISIBLE in background stats (z only -0.7..-1.8) and
            only reveal themselves when real energy arrives (z collapses to -6.5);
          - a MEDIAN over excited windows dilutes the signature away (caught 6/18);
          - LOCAL (±30 ch) normalisation fails completely on this array (caught 0/18):
            the weak zone is a ~50-channel contiguous cluster which normalises itself
            away. Global per-window z with an any-window vote catches 15/18 + 3
            same-zone borderline channels. Hence: global z, vote per excited window.
      * weak responders are ADVISORY by default (weak_in_mask=False): they carry good
        background and DO trigger on strong events (measured: all 18 were linked into
        the 15:09:33 chain) — statically masking them would discard those real picks.
        The picker's per-event SNR gate is the correct layer for them. weak_in_mask
        exists for the documented A/B experiment.

Window scheme: everything is computed on a scan_window_s (default 10 s) grid — long
windows dilute 10-30 s events. A window is EXCITED when the array-median log-energy
exceeds the scan median by max(excite_k_mad·1.4826·MAD, excite_min_decades); the rest
are BACKGROUND. This split is what keeps a passing ship or an earthquake from ever
making a channel "noisy".

Artifacts (one folder per scan, under <output_dir>/ChannelQuality/scan_<...>/):
    channel_quality.csv   one row per channel: all metrics + flags + mask + reason
    channel_mask.csv      minimal consumable (channel, dist_km, masked, reason) — hand-editable
    channel_quality.h5    per-window energy matrix E + metadata: re-threshold instantly
                          (run_v3.py mask) without touching the raw DAS
    quality_report.png    4-panel diagnostic figure
    scan_config.json      provenance (data dir, band, range, thresholds, time range)
"""

import os
import json
import glob
import numpy as np
from datetime import timedelta

from stalta import sta_lta_centred, trig_onechan_and_outindex     # vendored (v2 path)
from .qcio import read_das_qc

try:
    import h5py
except Exception:
    h5py = None

_EPS = 1e-15


# ====================================================================== scan ==========
def run_scan(cfg, progress=True):
    """Scan cfg.time_range, write artifacts, return (scan_dir, stats dict).

    Reads data through read_das_qc — identical loading + preprocessing to v3 detection
    (and to v2, when the raw data is clean).
    """
    q = cfg.quality
    win_s = float(q["scan_window_s"])
    scan_dir = cfg.scan_name(cfg.t_start, cfg.t_end)
    os.makedirs(scan_dir, exist_ok=True)

    chunks = []
    t = cfg.t_start
    while t < cfg.t_end:
        t2 = min(t + cfg.dt_chunk, cfg.t_end)
        chunks.append((t, t2))
        t = t2

    E_rows, wt_rows = [], []                       # per-window log10 RMS [nch], start time
    dead_raw_any, trig_total, dur_total, nch = None, None, 0.0, None

    for i, (a, b) in enumerate(chunks):
        if progress:
            print(f"[scan {i+1}/{len(chunks)}] {a} -> {b}")
        data, dead_raw = read_das_qc(a, b, cfg)
        nt = data.shape[0]
        if nch is None:
            nch = data.shape[1]
            trig_total = np.zeros(nch, dtype=int)
            dead_raw_any = np.zeros(nch, dtype=bool)
        dead_raw_any |= dead_raw

        # trigger rate at the configured detection parameters (report/corroborator)
        for k in range(nch):
            with np.errstate(all="ignore"):
                csl = sta_lta_centred(data[:, k], cfg.nsta, cfg.nlta)
            trig_total[k] += len(trig_onechan_and_outindex(
                np.nan_to_num(csl), cfg.threshold_level,
                cfg.minimal_event_interval, cfg.fs))
        dur_total += nt / cfg.fs

        nwin = int(nt // (win_s * cfg.fs))
        for w in range(nwin):
            d = data[int(w * win_s * cfg.fs): int((w + 1) * win_s * cfg.fs)]
            E_rows.append(np.log10(np.sqrt(np.mean(d ** 2, axis=0)) + _EPS))
            wt_rows.append((a + timedelta(seconds=w * win_s)).isoformat())

    E = np.asarray(E_rows)                          # [nwin, nch]
    if E.size == 0:
        raise RuntimeError(
            f"scan produced no windows: the time range ({dur_total:.0f} s of data) is "
            f"shorter than one scan window ({win_s:g} s). Extend time_range or lower "
            f"quality.scan_window_s.")
    trig_hr = trig_total * 3600.0 / max(dur_total, 1e-9)
    dist_km = (cfg.d1 + np.arange(nch) * cfg.dx) / 1000.0

    if h5py is not None:
        with h5py.File(os.path.join(scan_dir, "channel_quality.h5"), "w") as f:
            f.create_dataset("E_window", data=E.astype(np.float32), compression="gzip")
            f.create_dataset("window_start", data=np.array(wt_rows, dtype="S32"))
            f.create_dataset("dead_raw", data=dead_raw_any)
            f.create_dataset("trig_per_hr", data=trig_hr)
            f.create_dataset("dist_km", data=dist_km)
            f.attrs["scan_seconds"] = dur_total
            f.attrs["fs"] = cfg.fs

    stats = derive_mask(E, dead_raw_any, trig_hr, q)
    _write_artifacts(scan_dir, cfg, stats, dist_km, E)
    if progress:
        print_scan_summary(cfg, stats, nch)
        print(f">>> artifacts -> {scan_dir}")
    return scan_dir, stats


def print_scan_summary(cfg, stats, nch):
    q = cfg.quality
    m = stats["mask"]
    print(f">>> {stats['n_windows']} windows ({q['scan_window_s']:g} s): "
          f"{int(stats['hi'].sum())} excited, {int((~stats['hi']).sum())} background")
    weak_note = (f"weak {int(stats['weak'].sum())} "
                 f"(≥{stats.get('votes_needed', 1)} votes/"
                 f"{int(stats['hi'].sum())} excited; "
                 f"{'IN mask' if q.get('weak_in_mask') else 'advisory only'})")
    print(f">>> masked {int(m.sum())}/{nch} channels — "
          f"dead {int(stats['dead'].sum())} · noisy {int(stats['noisy'].sum())} · "
          f"{weak_note} · trigger-rate report-only")
    if not stats["hi"].any():
        print(">>> WARNING: no excited (high-energy) windows in the scan period — "
              "weak-responder detection UNMEASURED. Scan a period containing strong "
              "events to measure responsiveness.")
    for w in mask_guardrails(m, cfg):
        print(f">>> WARNING: {w}")


def mask_guardrails(mask, cfg):
    """Sanity warnings on a mask: excessive fraction; contiguous runs the tracker
    cannot bridge (run > neighbor_tol-1 splits every crossing wavefront)."""
    warns = []
    if mask is None or not len(mask):
        return warns
    frac = float(np.mean(mask))
    frac_alarm = float(cfg.quality.get("mask_frac_alarm", 0.15))
    if frac > frac_alarm:
        warns.append(f"mask covers {frac:.0%} of the array (> {frac_alarm:.0%}) — "
                     "check thresholds before trusting the catalog")
    ntol = int(cfg.tracking_params.get("neighbor_tol", 10))
    run, best, best_at = 0, 0, -1
    for k, v in enumerate(mask):
        run = run + 1 if v else 0
        if run > best:
            best, best_at = run, k - run + 1
    if best > ntol - 1:
        warns.append(f"contiguous masked run of {best} channels at ch{best_at}-"
                     f"{best_at+best-1} exceeds neighbor_tol-1={ntol-1}: wavefronts "
                     "crossing it WILL split (merge_time_tol may reunite them — "
                     "verify with run_v3.py compare)")
    return warns


# ====================================================================== flags =========
def derive_mask(E, dead_raw, trig_hr, q):
    """Pure function: window energy matrix + raw flags -> flags + mask.

    Used by scan and by re-threshold (run_v3.py mask) — no DAS access.
    """
    nwin, nch = E.shape

    # excitation split on the array-median energy
    M = np.median(E, axis=1)
    madM = np.median(np.abs(M - np.median(M)))
    thr = np.median(M) + max(float(q["excite_k_mad"]) * 1.4826 * madM,
                             float(q["excite_min_decades"]))
    hi = M >= thr

    # global robust z per window
    zg = np.empty_like(E)
    for w in range(nwin):
        med = np.median(E[w])
        mad = np.median(np.abs(E[w] - med)) + 1e-30
        zg[w] = (E[w] - med) / (1.4826 * mad)

    # background aggregation (static level) — excited windows excluded by construction
    bg = ~hi
    z_bg = np.median(zg[bg], axis=0) if bg.any() else np.full(nch, np.nan)

    # responsiveness: vote per excited window on the GLOBAL z (see module docstring).
    # The vote threshold scales with the number of excited windows (resp_min_frac):
    # calibrated on the full 2024-04-25 day scan — with 287 excited windows a single
    # vote is geometry/radiation-pattern chance (7 channels were flagged on 1/287),
    # while the genuine 8.6-9.0 km deaf belt fails 9-32% of excited windows.
    if hi.any():
        votes = (zg[hi] < float(q["resp_z"])).sum(axis=0)
        z_hi_min = np.min(zg[hi], axis=0)
    else:
        votes = np.zeros(nch, dtype=int)
        z_hi_min = np.full(nch, np.nan)
    votes_needed = max(int(q["resp_min_votes"]),
                       int(np.ceil(float(q.get("resp_min_frac", 0.03)) * int(hi.sum()))))

    # trigger-rate corroborator
    lr = np.log10(1.0 + trig_hr)
    med = np.median(lr)
    mad = np.median(np.abs(lr - med)) + 1e-30
    z_rate = (lr - med) / (1.4826 * mad)

    dead = dead_raw | ~np.isfinite(z_bg) | np.where(np.isfinite(z_bg),
                                                    z_bg <= float(q["z_dead"]), False)
    noisy = np.where(np.isfinite(z_bg), z_bg >= float(q["z_noisy"]), False) | (
        np.where(np.isfinite(z_bg), z_bg >= float(q["z_noisy_soft"]), False)
        & (z_rate >= float(q["z_rate"]))
        & (trig_hr >= float(q["rate_min_per_hr"])))
    weak = votes >= votes_needed if hi.any() else np.zeros(nch, bool)
    trig_happy = trig_hr > float(q["trig_hr_max"])          # report-only, never masks

    # local-excess z: background z minus an along-fiber running-median trend. REPORT
    # only — it separates a genuinely anomalous single channel from a regionally-noisy
    # environment (the 13.5-13.9 km tail is +2..+3 globally but ~0 vs its neighbours).
    # Not used for masks: a wide bad cluster would self-normalise (measured 0/18 on
    # the weak belt), the same trap that disqualified local z for responsiveness.
    from scipy.ndimage import median_filter
    win_t = int(q.get("local_trend_ch", 61)) | 1
    if np.isfinite(z_bg).all() and nch >= win_t:
        trend = median_filter(z_bg, size=win_t, mode="nearest")
        resid = z_bg - trend
        mad_r = np.median(np.abs(resid - np.median(resid))) + 1e-30
        z_local = (resid - np.median(resid)) / (1.4826 * mad_r)
    else:
        z_local = np.full(nch, np.nan)

    mask = dead | noisy
    if q.get("weak_in_mask", False):
        mask = mask | weak
    manual_bad = np.zeros(nch, bool)
    manual_good = np.zeros(nch, bool)
    for k in q.get("manual_bad", []):
        if 0 <= int(k) < nch:
            manual_bad[int(k)] = True
    for k in q.get("manual_good", []):
        if 0 <= int(k) < nch:
            manual_good[int(k)] = True
    mask = (mask | manual_bad) & ~manual_good

    weak_in_mask = bool(q.get("weak_in_mask", False))
    reasons = []
    for k in range(nch):
        r = [n for n, f in [("dead", dead), ("noisy", noisy),
                            ("weak_responder", weak if weak_in_mask else
                             np.zeros(nch, bool)),
                            ("manual", manual_bad)] if f[k]]
        if manual_good[k]:
            r = []
        reasons.append("+".join(r))

    return {"z_bg": z_bg, "z_local": z_local, "z_hi_min": z_hi_min, "votes": votes,
            "votes_needed": votes_needed, "trig_hr": trig_hr,
            "z_rate": z_rate, "hi": hi, "n_windows": nwin,
            "dead": dead, "noisy": noisy, "weak": weak, "trig_happy": trig_happy,
            "manual_bad": manual_bad, "manual_good": manual_good,
            "mask": mask, "reasons": reasons}


# ====================================================================== artifacts =====
def _write_artifacts(scan_dir, cfg, stats, dist_km, E):
    """Full artifact set for a fresh scan (tables + provenance + figure)."""
    _write_tables(scan_dir, stats, dist_km)
    prov = {
        "das_data_dir": cfg.das_data_dir, "fs": cfg.fs, "dx": cfg.dx,
        "freq_band": cfg.freq_band, "spatial_range": cfg.spatial_range,
        "time_range": [cfg.t_start.isoformat(), cfg.t_end.isoformat()],
        "sta_len": cfg.sta_len, "lta_len": cfg.lta_len,
        "threshold_level": cfg.threshold_level, "quality": cfg.quality,
        "n_channels": int(len(dist_km)), "n_windows": int(stats["n_windows"]),
        "n_excited_windows": int(stats["hi"].sum()),
        "n_masked": int(stats["mask"].sum()),
    }
    with open(os.path.join(scan_dir, "scan_config.json"), "w") as f:
        json.dump(prov, f, indent=2, ensure_ascii=False)
    try:
        _report_figure(os.path.join(scan_dir, "quality_report.png"),
                       cfg.quality, f"{cfg.t_start} → {cfg.t_end}",
                       stats, dist_km, E)
    except Exception as e:                            # the figure is best-effort
        print(f"[scan] report figure failed: {e}")


def rewrite_artifacts(scan_dir, quality, stats, dist_km, E):
    """Refresh tables + figure after a re-threshold, WITHOUT touching the scan's
    recorded data provenance (only the quality thresholds + mask count update)."""
    _write_tables(scan_dir, stats, dist_km)
    prov_path = os.path.join(scan_dir, "scan_config.json")
    prov, time_note = {}, "?"
    if os.path.exists(prov_path):
        with open(prov_path) as f:
            prov = json.load(f)
        tr = prov.get("time_range", ["?", "?"])
        time_note = f"{tr[0]} → {tr[1]}"
    prov["quality"] = quality
    prov["n_masked"] = int(stats["mask"].sum())
    with open(prov_path, "w") as f:
        json.dump(prov, f, indent=2, ensure_ascii=False)
    try:
        _report_figure(os.path.join(scan_dir, "quality_report.png"),
                       quality, time_note, stats, dist_km, E)
    except Exception as e:
        print(f"[mask] report figure failed: {e}")


def _write_tables(scan_dir, stats, dist_km):
    import pandas as pd
    nch = len(dist_km)
    df = pd.DataFrame({
        "channel": np.arange(nch),
        "dist_km": np.round(dist_km, 4),
        "z_bg": np.round(stats["z_bg"], 3),
        "z_local": np.round(stats["z_local"], 3),
        "z_hi_min": np.round(stats["z_hi_min"], 3),
        "resp_votes": stats["votes"],
        "trig_per_hr": np.round(stats["trig_hr"], 1),
        "z_rate": np.round(stats["z_rate"], 2),
        "flag_dead": stats["dead"], "flag_noisy": stats["noisy"],
        "flag_weak_responder": stats["weak"],
        "flag_trigger_happy": stats["trig_happy"],
        "manual_bad": stats["manual_bad"], "manual_good": stats["manual_good"],
        "masked": stats["mask"], "reason": stats["reasons"],
    })
    df.to_csv(os.path.join(scan_dir, "channel_quality.csv"), index=False)
    df[["channel", "dist_km", "masked", "reason"]].to_csv(
        os.path.join(scan_dir, "channel_mask.csv"), index=False)


def _report_figure(path, q, time_note, stats, dist_km, E):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    hi = stats["hi"]
    nwin, nch = E.shape
    # z matrix for display (global per-window z, same as the flags use)
    zg = np.empty_like(E)
    for w in range(nwin):
        med = np.median(E[w])
        mad = np.median(np.abs(E[w] - med)) + 1e-30
        zg[w] = (E[w] - med) / (1.4826 * mad)

    fig, axes = plt.subplots(1, 4, figsize=(17, 8), sharey=True,
                             gridspec_kw={"width_ratios": [2.2, 1, 1, 1]})
    im = axes[0].imshow(zg.T, aspect="auto", cmap="coolwarm", vmin=-5, vmax=5,
                        extent=[0, nwin, dist_km.max(), dist_km.min()],
                        interpolation="nearest")
    for w in np.where(hi)[0]:
        axes[0].axvline(w + 0.5, color="k", lw=0.5, alpha=0.30)
    axes[0].set_xlabel(f"scan window # ({q['scan_window_s']:g} s each; "
                       f"ticks = {int(hi.sum())} excited)")
    axes[0].set_ylabel("distance along cable (km)")
    axes[0].set_title("robust z of log10(RMS) per window")
    fig.colorbar(im, ax=axes[0], fraction=0.03)

    axes[1].plot(stats["z_bg"], dist_km, lw=0.8)
    axes[1].axvline(q["z_dead"], color="k", ls="--", lw=0.8)
    axes[1].axvline(q["z_noisy"], color="r", ls="--", lw=0.8)
    axes[1].set_xlabel("background z (median)")
    axes[1].set_title("static level\n(dead | noisy)")
    axes[1].grid(alpha=0.3)

    axes[2].plot(stats["z_hi_min"], dist_km, lw=0.8, color="purple")
    axes[2].axvline(q["resp_z"], color="r", ls="--", lw=0.8)
    axes[2].set_xlabel("min z in excited windows")
    axes[2].set_title("responsiveness\n(weak responder"
                      + ("" if q.get("weak_in_mask") else ", advisory") + ")")
    axes[2].grid(alpha=0.3)

    axes[3].plot(stats["trig_hr"], dist_km, lw=0.8, color="g")
    axes[3].axvline(q["trig_hr_max"], color="r", ls="--", lw=0.8)
    axes[3].set_xlabel("triggers / hour")
    axes[3].set_title("trigger rate\n(report-only)")
    axes[3].grid(alpha=0.3)

    for ax in axes[1:]:
        for k in np.where(stats["mask"])[0]:
            ax.axhspan(dist_km[k] - 0.006, dist_km[k] + 0.006,
                       color="orange", alpha=0.25, linewidth=0)
        for k in np.where(stats["weak"] & ~stats["mask"])[0]:
            ax.axhspan(dist_km[k] - 0.006, dist_km[k] + 0.006,
                       color="purple", alpha=0.12, linewidth=0)
    fig.suptitle(f"channel quality scan · {time_note} · "
                 f"{int(stats['mask'].sum())}/{nch} masked (orange) · "
                 f"{int((stats['weak'] & ~stats['mask']).sum())} weak advisory (purple)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ====================================================================== consumption ===
def rethreshold(scan_dir, quality):
    """Recompute flags + mask from a stored scan (no DAS read)."""
    if h5py is None:
        raise RuntimeError("h5py required")
    with h5py.File(os.path.join(scan_dir, "channel_quality.h5"), "r") as f:
        E = f["E_window"][()].astype(float)
        dead_raw = f["dead_raw"][()]
        trig_hr = f["trig_per_hr"][()]
        dist_km = f["dist_km"][()]
    stats = derive_mask(E, dead_raw, trig_hr, quality)
    return stats, dist_km, E


def load_mask(path, nch=None):
    """Read a channel_mask.csv -> (bool mask [nch], reasons list). Tolerates hand edits."""
    import pandas as pd
    df = pd.read_csv(path)
    n = int(df["channel"].max()) + 1 if nch is None else int(nch)
    mask = np.zeros(n, bool)
    reasons = [""] * n
    for _, r in df.iterrows():
        k = int(r["channel"])
        if 0 <= k < n:
            m = r["masked"]
            if isinstance(m, str):
                mask[k] = m.strip().lower() in ("true", "1", "yes")
            elif pd.isna(m):                 # blanked cell in a hand-edited csv = unmask
                mask[k] = False
            else:
                mask[k] = bool(m)
            reasons[k] = r["reason"] if isinstance(r.get("reason", ""), str) else ""
    return mask, reasons


def find_latest_mask(cfg):
    """Auto-discover the most recent compatible scan's channel_mask.csv (or None).

    Compatible = same data dir, fs, band and spatial range (via scan_config.json)."""
    pattern = os.path.join(cfg.quality_root, "scan_*", "channel_mask.csv")
    cands = sorted(glob.glob(pattern), key=os.path.getmtime, reverse=True)
    for c in cands:
        try:
            with open(os.path.join(os.path.dirname(c), "scan_config.json")) as f:
                prov = json.load(f)
            if (prov.get("das_data_dir") == cfg.das_data_dir
                    and float(prov.get("fs", -1)) == float(cfg.fs)
                    and float(prov.get("dx", -1)) == float(cfg.dx)
                    and list(map(float, prov.get("freq_band", []))) == [cfg.f1, cfg.f2]
                    and list(map(float, prov.get("spatial_range", []))) == [cfg.d1, cfg.d2]):
                return c
        except Exception:
            continue
    return None


def resolve_mask(cfg, no_mask=False, require=False):
    """Attach the channel mask to cfg (cfg.chan_mask / chan_mask_path / chan_mask_reasons).

    Order: --no-mask -> None;  quality.mask_file -> that file;  else auto-discover.
    An explicitly-given file whose sibling scan_config.json disagrees with the live
    geometry triggers a loud warning (indices would mean different channels).
    """
    if no_mask:
        cfg.chan_mask, cfg.chan_mask_path = None, ""
        return cfg
    nch = int(cfg.d2 / cfg.dx) - int(cfg.d1 / cfg.dx)
    path = cfg.quality.get("mask_file") or ""
    if not path:
        path = find_latest_mask(cfg) or ""
    else:
        prov_path = os.path.join(os.path.dirname(os.path.abspath(path)),
                                 "scan_config.json")
        if os.path.exists(prov_path):
            try:
                with open(prov_path) as f:
                    prov = json.load(f)
                if (list(map(float, prov.get("spatial_range", []))) != [cfg.d1, cfg.d2]
                        or float(prov.get("fs", cfg.fs)) != float(cfg.fs)
                        or float(prov.get("dx", cfg.dx)) != float(cfg.dx)
                        or list(map(float, prov.get("freq_band",
                                                    [cfg.f1, cfg.f2]))) != [cfg.f1, cfg.f2]):
                    print(f"[v3] WARNING: mask {path} was scanned with a different "
                          f"geometry (fs/dx/band/spatial_range) — channel indices may "
                          f"not line up with this configuration.")
            except Exception:
                pass
    if not path:
        if require:
            raise RuntimeError(
                "no channel mask found — run `python run_v3.py scan` first "
                "(or pass --mask <channel_mask.csv> / --no-mask).")
        cfg.chan_mask, cfg.chan_mask_path = None, ""
        return cfg
    mask, reasons = load_mask(path, nch)
    try:                                      # row-count sanity for hand-supplied files
        import pandas as pd
        n_rows = int(pd.read_csv(path)["channel"].max()) + 1
        if n_rows != nch:
            print(f"[v3] WARNING: mask {path} covers {n_rows} channels but the current "
                  f"geometry has {nch} — missing channels default to unmasked; check "
                  f"that the mask matches this fs/dx/spatial_range.")
    except Exception:
        pass
    cfg.chan_mask = mask
    cfg.chan_mask_path = os.path.abspath(path)
    cfg.chan_mask_reasons = reasons
    return cfg
