"""
Unified configuration for das_wavefront_picker.

ONE config object drives the whole app (detection + interactive picking + export),
so switching between the bundled 80 Hz / 16 m demo and 400 Hz / 8 m production data is
just a matter of changing values in `config.json` (or passing --flags). No environment
variables are used; every path/parameter is explicit and visible.

Layout assumed:
    das_wavefront_picker/
        vendor/functions/         # vendored STA/LTA kernel + DAS reader (self-contained)
        dasapp/config.py          # <- this file
        config.json               # user-editable defaults
"""

import os
import sys
import glob
import json
import argparse
from datetime import datetime, timedelta

# ---- locate the app root and put the vendored kernel on sys.path --------------------
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
APP_ROOT = os.path.dirname(_THIS_DIR)
VENDOR_FUNCTIONS = os.path.join(APP_ROOT, "vendor", "functions")
if VENDOR_FUNCTIONS not in sys.path:
    sys.path.insert(0, VENDOR_FUNCTIONS)


def _find_demo_data_dir():
    """Best-effort locate the demo DAS sample so `run.py demo` / `selftest` work without
    editing paths. Checks, in order:
      1. a `demo_data/` folder dropped INSIDE this app (the portable way — get the sample
         h5 from matched_filter_release/IN_DAS_data and put them in das_wavefront_picker/demo_data/);
      2. the matched_filter_release location in the original dev tree.
    Returns None if no folder with .h5 files is found (demo/selftest then print a clear hint)."""
    candidates = [
        os.path.join(APP_ROOT, "demo_data"),
        os.path.normpath(os.path.join(
            APP_ROOT, "..", "..", "PRE_DISSCUSSION", "Kefa-DAS-GRL",
            "matched_filter_release", "stalta_multichan_v1", "IN_DAS_data")),
    ]
    for c in candidates:
        if os.path.isdir(c) and glob.glob(os.path.join(c, "**", "*.h5"), recursive=True):
            return c + os.sep
    return None


DEMO_DATA_DIR = _find_demo_data_dir()

# ---- built-in defaults (bundled 80 Hz / 16 m demo, 2024-04-25 15:00-15:30) ----------
DEFAULTS = {
    # data + preprocessing
    "das_data_dir": DEMO_DATA_DIR or "",
    "fs": 80.0,
    "dx": 16 * 1.02095199,          # channel spacing (m)
    "freq_band": [5.0, 20.0],       # bandpass corners (Hz)
    "spatial_range": [8000.0, 14000.0],   # d1, d2 (m) region of interest

    # STA/LTA
    "sta_len": 0.3,                 # short-term window (s)
    "lta_len": 3.0,                 # long-term window (s)

    # detection
    "method": "B",                  # "A" (simple count) or "B" (wavefront tracking; seeds picks)
    "threshold_level": 35.0,        # single-channel STA/LTA trigger threshold (MAD-normalised)
    "minimal_event_interval": 1.2,  # debounce interval on one channel (s)
    "min_channels": 50,             # min channels for an event (trig_NO / threshold_channels_per_event)
    "tracking_params": {"neighbor_tol": 10, "time_tol": 0.05, "slope_tol": 0.05},  # method B
    "merge_params": {"merge_time_tol": 1.0},                                        # method B
    "time_range": ["2024-04-25T15:00:00", "2024-04-25T15:30:00"],
    "dt_chunk_min": 10.0,           # processing chunk length (minutes)

    # output
    "output_dir": "Output",         # relative to app root unless absolute

    # interactive picking
    "labels": ["P", "S", "vehicle", "explosive", "noise", "unknown"],  # controlled vocabulary
    "default_label": "unknown",     # label a freshly-loaded wavefront starts with
    "pad_pre_s": 2.0,               # raw-DAS window pad before the candidate reference time (s)
    "pad_post_s": 3.0,              # ... after event end (s)
    "n_anchors": 16,                # number of movable anchor control points along the array
    "n_wiggle": 30,                 # initial traces in the wiggle background (live-tunable with , / .)
    "edge_margin_km": 0.0,          # channels within this margin of the array ends -> low-confidence
    "wavefront_alpha": 0.6,         # opacity of the green wavefront line (lower = see more background)
}


class Config:
    """Loaded from DEFAULTS <- config.json <- CLI flags (later overrides earlier)."""

    def __init__(self, overrides=None, config_path=None):
        cfg = dict(DEFAULTS)

        # merge config.json (next to app root) if present
        self.config_path = config_path or os.path.join(APP_ROOT, "config.json")
        if os.path.exists(self.config_path):
            with open(self.config_path) as f:
                disk = json.load(f)
            cfg.update({k: v for k, v in disk.items() if v is not None})

        # merge explicit overrides (CLI)
        if overrides:
            cfg.update({k: v for k, v in overrides.items() if v is not None})

        # store
        for k, v in cfg.items():
            setattr(self, k, v)

        # resolve paths
        self.das_data_dir = self._resolve(self.das_data_dir) if self.das_data_dir else ""
        self.output_dir = self._resolve(self.output_dir)

        # parse time range
        self.t_start, self.t_end = (self._parse_dt(x) for x in self.time_range)
        self.dt_chunk = timedelta(minutes=float(self.dt_chunk_min))

    # ---- derived ----
    @property
    def nsta(self):
        return int(self.sta_len * self.fs)          # truncate (matches the release loops)

    @property
    def nlta(self):
        return int(self.lta_len * self.fs)          # truncate (matches the release loops)

    @property
    def d1(self):
        return float(self.spatial_range[0])

    @property
    def d2(self):
        return float(self.spatial_range[1])

    @property
    def f1(self):
        return float(self.freq_band[0])

    @property
    def f2(self):
        return float(self.freq_band[1])

    @property
    def exp_name(self):
        # include times, not just dates, so two intra-day windows don't collide into one folder
        t = f"{self.t_start.strftime('%Y%m%d-%H%M%S')}_{self.t_end.strftime('%Y%m%d-%H%M%S')}"
        return (f"method{self.method}_thr{self.threshold_level:g}_minCh{self.min_channels}"
                f"_sta{self.sta_len:g}_lta{self.lta_len:g}_{t}")

    @property
    def exp_dir(self):
        return os.path.join(self.output_dir, self.exp_name)

    @property
    def catalog_path(self):
        return os.path.join(self.exp_dir, "catalog.csv")

    @property
    def wavefront_dir(self):
        return os.path.join(self.exp_dir, "Wavefronts")

    @property
    def picks_dir(self):
        return os.path.join(self.exp_dir, "Picks")

    # ---- helpers ----
    @staticmethod
    def _resolve(p):
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(APP_ROOT, p))

    @staticmethod
    def _parse_dt(x):
        if isinstance(x, datetime):
            return x
        return datetime.fromisoformat(str(x))

    def to_dict(self):
        keys = ["das_data_dir", "fs", "dx", "freq_band", "spatial_range", "sta_len",
                "lta_len", "method", "threshold_level", "minimal_event_interval",
                "min_channels", "tracking_params", "merge_params", "time_range",
                "dt_chunk_min", "output_dir", "labels", "default_label",
                "pad_pre_s", "pad_post_s", "n_anchors", "n_wiggle", "edge_margin_km",
                "wavefront_alpha"]
        return {k: getattr(self, k) for k in keys}

    def save_snapshot(self, path=None):
        path = path or os.path.join(self.exp_dir, "config_used.json")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)
        return path

    def summary(self):
        return (f"data={self.das_data_dir}\n"
                f"fs={self.fs} dx={self.dx:.3f}m band={self.freq_band}Hz "
                f"range={self.spatial_range}m\n"
                f"STA/LTA={self.sta_len}/{self.lta_len}s thr={self.threshold_level} "
                f"minCh={self.min_channels} method={self.method}\n"
                f"time={self.t_start} .. {self.t_end}\n"
                f"out={self.exp_dir}")


def add_common_cli(ap):
    """Attach the shared --flags that override config.json for any subcommand."""
    ap.add_argument("--config", help="path to a config.json (default: app root config.json)")
    ap.add_argument("--data", dest="das_data_dir", help="raw DAS .h5 directory")
    ap.add_argument("--fs", type=float)
    ap.add_argument("--dx", type=float)
    ap.add_argument("--band", nargs=2, type=float, metavar=("F1", "F2"), dest="freq_band")
    ap.add_argument("--range", nargs=2, type=float, metavar=("D1", "D2"), dest="spatial_range")
    ap.add_argument("--sta", type=float, dest="sta_len")
    ap.add_argument("--lta", type=float, dest="lta_len")
    ap.add_argument("--method", choices=["A", "B"])
    ap.add_argument("--threshold", type=float, dest="threshold_level")
    ap.add_argument("--min-channels", type=int, dest="min_channels")
    ap.add_argument("--start", dest="t_start_str", help="ISO start, e.g. 2024-04-25T15:00:00")
    ap.add_argument("--end", dest="t_end_str", help="ISO end")
    ap.add_argument("--out", dest="output_dir")


def config_from_args(args):
    """Build a Config from parsed argparse args (only non-None flags override)."""
    ov = {}
    for k in ["das_data_dir", "fs", "dx", "freq_band", "spatial_range", "sta_len",
              "lta_len", "method", "threshold_level", "min_channels", "output_dir"]:
        v = getattr(args, k, None)
        if v is not None:
            ov[k] = v
    cfg = Config(overrides=ov, config_path=getattr(args, "config", None))
    # time-range overrides: --start and --end apply independently (either one alone is honoured,
    # falling back to the config value for the side not given)
    ts = getattr(args, "t_start_str", None)
    te = getattr(args, "t_end_str", None)
    if ts or te:
        a = ts or cfg.time_range[0]
        b = te or cfg.time_range[1]
        cfg.time_range = [a, b]
        cfg.t_start = Config._parse_dt(a)
        cfg.t_end = Config._parse_dt(b)
    return cfg
