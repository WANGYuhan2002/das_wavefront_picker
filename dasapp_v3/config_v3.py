"""
v3 configuration: the v2 Config plus a `quality` section (channel-quality scan + mask).

Design rules (see V3_CHANGES.md):
  * dasapp/config.py is NOT modified — ConfigV3 subclasses it.
  * All quality parameters live in ONE nested dict so config.json stays tidy:
        "quality": { "z_hot": 3.0, ... }        (any subset; the rest use defaults)
  * v3 detection output NEVER collides with v2 output: the experiment folder name is
    prefixed "v3_" and carries a short hash of the channel mask actually used.

Channel indexing convention: everywhere in the app (and in every v3 artifact) a
"channel" is the index INSIDE the spatial slice int(d1/dx):int(d2/dx) — the same
convention the detector, picker and wavefront H5 files already use. Artifacts also
store dist_km so files remain interpretable if spatial_range changes.
"""

import os
import hashlib

from dasapp.config import Config, add_common_cli, DEMO_DATA_DIR  # noqa: F401 (re-export)

# ---- quality defaults (all overridable via config.json "quality": {...} or CLI) ------
# Every threshold below was calibrated on the bundled demo; the measurements and the
# reasoning live in V3_CHANGES.md. z values are robust z-scores of log10(RMS) across
# channels (MAD-normalised), computed per scan window.
QUALITY_DEFAULTS = {
    # ---- static scan (run_v3.py scan) ----
    "scan_window_s": 10.0,      # window grid (s); long windows dilute 10-30 s events
    "excite_k_mad": 3.0,        # a window is EXCITED when the array-median log-energy
    "excite_min_decades": 0.5,  #   exceeds the scan median by max(k_mad·1.4826·MAD, decades)
    "z_dead": -4.0,             # DEAD: background z at/below this (plus raw-domain screen)
    "z_noisy": 3.0,             # NOISY: background z at/above this...
    "z_noisy_soft": 2.5,        # ...or above this WITH trigger-rate corroboration:
    "z_rate": 3.5,              #    robust z of log10(1+trig/hr) at/above this...
    "rate_min_per_hr": 20.0,    #    ...and an absolute rate of at least this
    "resp_z": -3.0,             # WEAK RESPONDER: global z below this in an excited window...
    "resp_min_votes": 1,        # ...in at least max(this, ceil(resp_min_frac·n_excited))
    "resp_min_frac": 0.03,      #    windows. Fraction rule calibrated on the full-day scan:
                                #    287 excited windows -> need >=9 votes; the genuine deaf
                                #    belt fails 9-32%, single-vote flukes are excluded.
    "local_trend_ch": 61,       # running-median window (channels) for the REPORT-ONLY
                                #    local-excess column z_local (regional vs local noise)
    "weak_in_mask": False,      # weak responders are ADVISORY by default: they see strong
                                #   events (measured) and the picker's per-event gate is the
                                #   right layer; True = the documented A/B experiment
    "trig_hr_max": 60.0,        # trigger-rate REPORT threshold (never masks by itself:
                                #   on active data a busy channel is a GOOD channel)
    "mask_frac_alarm": 0.15,    # warn when the mask covers more than this fraction
    "manual_bad": [],           # channel indices always masked (user override)
    "manual_good": [],          # channel indices never masked (wins over every flag)
    # ---- detection consumption ----
    "scale_min_channels": True, # min_channels -> round(min_channels · n_unmasked/n_total)
    # ---- event-local gate (picker v3) ----
    "event_snr_gate": 0.5,      # unlinked channels need local CF peak > gate·threshold_level
                                #   (same convention as v2 picker._seed_auto for method-A)
    "event_snr_win_s": 0.5,     # ± window (s) around the seed arrival for the local CF peak
    "prune_k": 3.0,             # linked channels: demote if |arrival − local median| >
                                #   max(prune_k·1.4826·MAD, 3 samples) over prune_win_ch
    "prune_win_ch": 15,         # moving window (channels) for the mislink prune (0 = off)
    "snr_noise_s": [-2.5, -0.5],  # windows (s, rel. to pick) for the stored snr_db
    "snr_sig_s": [0.0, 0.5],      #   attribute (diagnostic only — never a gate)
    # ---- mask consumption ----
    "mask_file": "",            # explicit channel_mask.csv path; "" = auto-discover the most
                                #   recent compatible scan under <output_dir>/ChannelQuality/
}


class ConfigV3(Config):
    """v2 Config + quality section. `quality` merges DEFAULTS <- config.json <- overrides."""

    def __init__(self, overrides=None, config_path=None, quality_overrides=None):
        super().__init__(overrides=overrides, config_path=config_path)
        q = dict(QUALITY_DEFAULTS)
        disk_q = getattr(self, "quality", None)      # from config.json, if present
        if isinstance(disk_q, dict):
            q.update({k: v for k, v in disk_q.items() if v is not None})
        if quality_overrides:
            q.update({k: v for k, v in quality_overrides.items() if v is not None})
        self.quality = q

        # populated by run_v3 once a mask file is resolved/loaded:
        self.chan_mask = None            # bool array [nch] or None (no mask)
        self.chan_mask_path = ""         # provenance
        self.chan_mask_reasons = None    # list[str] per channel ("" if unmasked)

    # ---- mask-aware experiment naming (v3 outputs never collide with v2) -------------
    @property
    def mask_tag(self):
        if self.chan_mask is None or not self.chan_mask.any():
            return "nomask"
        h = hashlib.md5(self.chan_mask.tobytes()).hexdigest()[:6]
        return f"m{int(self.chan_mask.sum())}-{h}"

    @property
    def exp_name(self):
        return f"v3_{super().exp_name}_{self.mask_tag}"

    # ---- quality artifact locations --------------------------------------------------
    @property
    def quality_root(self):
        return os.path.join(self.output_dir, "ChannelQuality")

    def scan_name(self, t0, t1):
        """One scan = one folder, keyed by everything the statistics depend on."""
        tag = (f"{t0.strftime('%Y%m%d-%H%M%S')}_{t1.strftime('%Y%m%d-%H%M%S')}"
               f"_band{self.f1:g}-{self.f2:g}_range{self.d1:g}-{self.d2:g}"
               f"_win{self.quality['scan_window_s']:g}s")
        return os.path.join(self.quality_root, f"scan_{tag}")

    def to_dict(self):
        d = super().to_dict()
        d["quality"] = self.quality
        d["chan_mask_path"] = self.chan_mask_path
        if self.chan_mask is not None:
            d["chan_mask_n"] = int(self.chan_mask.sum())
            d["chan_mask_channels"] = [int(k) for k in self.chan_mask.nonzero()[0]]
        return d


def add_v3_cli(ap):
    """v3-specific flags, added on top of the shared v2 add_common_cli set."""
    ap.add_argument("--mask", dest="mask_file",
                    help="channel_mask.csv to use (default: auto-discover latest scan)")
    ap.add_argument("--no-mask", action="store_true",
                    help="run without any channel mask (pure v2 behaviour)")
    ap.add_argument("--scan-window", type=float, dest="scan_window_s",
                    help="quality scan window length (s)")
    ap.add_argument("--z-dead", type=float, dest="z_dead")
    ap.add_argument("--z-noisy", type=float, dest="z_noisy")
    ap.add_argument("--resp-z", type=float, dest="resp_z")
    ap.add_argument("--weak-in-mask", action="store_true", dest="weak_in_mask",
                    default=None, help="ALSO statically mask weak responders (A/B experiment)")
    ap.add_argument("--snr-gate", type=float, dest="event_snr_gate",
                    help="picker event-local conf gate, in units of threshold_level")


def config_v3_from_args(args):
    """Build a ConfigV3 from argparse args (v2 flags + v3 quality flags)."""
    ov = {}
    for k in ["das_data_dir", "fs", "dx", "freq_band", "spatial_range", "sta_len",
              "lta_len", "method", "threshold_level", "min_channels", "output_dir"]:
        v = getattr(args, k, None)
        if v is not None:
            ov[k] = v
    qov = {}
    for k in ["mask_file", "scan_window_s", "z_dead", "z_noisy", "resp_z",
              "weak_in_mask", "event_snr_gate"]:
        v = getattr(args, k, None)
        if v is not None:
            qov[k] = v
    cfg = ConfigV3(overrides=ov, config_path=getattr(args, "config", None),
                   quality_overrides=qov)
    ts = getattr(args, "t_start_str", None)
    te = getattr(args, "t_end_str", None)
    if ts or te:
        a = ts or cfg.time_range[0]
        b = te or cfg.time_range[1]
        cfg.time_range = [a, b]
        cfg.t_start = Config._parse_dt(a)
        cfg.t_end = Config._parse_dt(b)
    return cfg
