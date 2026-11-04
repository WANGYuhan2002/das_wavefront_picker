"""
Interactive DETECTION-parameter tuner (separate from picking).

Pick a representative window from an existing catalog, then drag sliders for all the Method-B
detection parameters and watch the STA/LTA + wavefront-tracking result update LIVE on that one
window. When it looks right, "Save to config" writes the parameters back to config.json — then run
the full batch `python run.py detect` with them.

Why a separate mode (not mixed into picking): detection decides *which events exist* — a global,
reproducible, batch decision recorded in config_used.json. Tuning it live inside the per-event
picker would make the catalog inconsistent and blur the detect/pick separation. Here we only ever
re-run on ONE preview window (cheap) and emit parameters; the catalog stays a single batch product.

Recompute cost is kept low: the STA/LTA characteristic function (csl) is cached and only recomputed
when STA/LTA change; the chain-tracking is re-run for the other parameters; min_channels is a pure
post-filter (no recompute). A single preview window is a *starting point* — check a few (strong /
weak / noisy) events with Prev/Next before committing, and confirm with one real detect run.
"""

import os
import json
import numpy as np
import matplotlib
from datetime import timedelta
from scipy.signal import hilbert

from stalta import (                                    # vendored
    sta_lta_centred_multichan,
    trig_allchan_and_collapsed_index_and_save_wavefonts,
)
from .dasio import read_window
from .picker import _FULLWIDTH, draw_shortcut_overlay    # CJK-IME key folding + the h-overlay

_TUNE_HELP_ROWS = [
    ("Background", "1 2 3 4", "raw amplitude · STA-LTA CF · envelope · wiggle traces"),
    ("Intensity",  "- / =",   "dim / brighten active bg (colour range · wiggle gain)"),
    ("Wiggle #",   ", / .",   "wiggle trace density: fewer / more"),
    ("Sliders",    "drag",    "the 9 detection params (meanings listed on the right)"),
    ("Buttons",    "click",   "Save to config · Reset params · Prev / Next event"),
    ("Help",       "h / q",   "toggle this help  ·  quit"),
]

_LEAD = 6.0            # seconds of lead loaded for CF (LTA) warm-up (covers lta_len up to 6 s)
_PAD_PRE = 6.0         # display window before the event
_PAD_POST = 10.0       # ... and after (a bit of context around the event)

# key, slider label, min, max, step, recompute-level, one-line meaning (English, renders anywhere)
PARAMS = [
    ("threshold_level",       "threshold",      5.0, 100.0, 1.0,   "track",
     "single-channel STA/LTA trigger level; higher -> fewer false triggers"),
    ("sta_len",               "STA (s)",        0.05, 1.0,  0.05,  "cf",
     "short-term window; sensitivity to impulsive onsets"),
    ("lta_len",               "LTA (s)",        0.5,  6.0,  0.1,   "cf",
     "long-term window; the background-noise level"),
    ("minimal_event_interval","min interval(s)",0.2,  3.0,  0.1,   "track",
     "debounce: same-channel triggers within this merge into one"),
    ("min_channels",          "min channels",   4,    800,  1,     "filter",
     "min channels a wavefront needs to count as an event (post-filter)"),
    ("neighbor_tol",          "neighbor tol",   1,    30,   1,     "track",
     "how many dead channels the tracking may skip"),
    ("time_tol",              "time tol (s)",   0.005,0.3,  0.005, "track",
     "max arrival-time gap between adjacent channels (max moveout)"),
    ("slope_tol",             "slope tol",      0.005,0.3,  0.005, "track",
     "wavefront slope-change tolerance (smoothness constraint)"),
    ("merge_time_tol",        "merge tol (s)",  0.1,  3.0,  0.1,   "track",
     "merge two wavefronts if their mean times differ by less than this"),
]


def _cfg_value(cfg, key):
    """Read a flat parameter's current value out of the nested config."""
    if key in ("neighbor_tol", "time_tol", "slope_tol"):
        return cfg.tracking_params.get(key)
    if key == "merge_time_tol":
        return cfg.merge_params.get("merge_time_tol")
    return getattr(cfg, key)


class Tuner:
    def __init__(self, cfg, cands, start_idx=0, interactive=True):
        self.cfg = cfg
        self.cands = cands
        self.interactive = interactive
        self.vals = {k: float(_cfg_value(cfg, k)) for (k, *_rest) in PARAMS}
        self.idx = int(np.clip(start_idx, 0, len(cands) - 1))
        self.fs = cfg.fs
        self.dx = cfg.dx
        self.d1 = cfg.d1
        self.csl = None
        self.wavefronts = []
        self.status = ""

        if not interactive:
            matplotlib.use("Agg")
        else:
            for _bk in ("MacOSX", "QtAgg", "TkAgg"):
                try:
                    matplotlib.use(_bk)
                    break
                except Exception:
                    continue
        import matplotlib.pyplot as plt
        for _k in list(plt.rcParams):
            if _k.startswith("keymap."):
                plt.rcParams[_k] = []
        self.plt = plt

        self.bg = ["2"]                                # 1 amplitude(raw) / 2 STA-LTA CF / 3 env / 4 wiggle
        self.n_wig = int(cfg.n_wiggle)                  # wiggle trace count (, / . adjust)
        self.gain = {"1": 1.0, "2": 1.0, "3": 1.0, "4": 1.0}   # per-bg intensity (- / = adjust)
        self.help_visible = [False]
        self.fig = plt.figure(figsize=(16, 10))
        self.ax = self.fig.add_axes([0.055, 0.45, 0.90, 0.47])   # waterfall: full width on top
        self._build_legend()
        self._build_sliders()
        self._build_buttons()
        if interactive:
            self.fig.canvas.mpl_connect("key_press_event", self._on_key)

        self._load_window()
        self._recompute("cf")

    def _on_key(self, ev):
        k = _FULLWIDTH.get(ev.key, ev.key)             # tolerate a CJK IME's full-width keys
        if k in ("1", "2", "3", "4"):                  # background: raw / CF / envelope / wiggle
            self.bg[0] = k
            self._draw()
        elif k in ("=", "+", "-"):                     # intensity of the active background
            f = 1.3 if k in ("=", "+") else 1 / 1.3
            self.gain[self.bg[0]] = float(np.clip(self.gain[self.bg[0]] * f, 0.1, 12.0))
            self._draw()
        elif k in (",", "."):                          # wiggle trace density
            f = 1.4 if k == "." else 1 / 1.4
            self.n_wig = int(np.clip(round(self.n_wig * f), 4, self.nch))
            if self.bg[0] == "4":
                self._draw()
        elif k == "h":                                 # toggle the shortcut overlay
            self.help_visible[0] = not self.help_visible[0]
            self._draw()
        elif k == "q":
            self.plt.close(self.fig)

    # ---- static meanings panel (bottom-right, above the buttons) ---------------------
    def _build_legend(self):
        lax = self.fig.add_axes([0.55, 0.205, 0.44, 0.205])
        lax.axis("off")
        lines = ["DETECTION PARAMETERS   ·   drag a slider → re-runs on this window", ""]
        for (k, label, *_r, meaning) in PARAMS:
            lines.append(f"{label}:  {meaning}")
        lax.text(0.0, 1.0, "\n".join(lines), transform=lax.transAxes, va="top", ha="left",
                 fontsize=8, family="DejaVu Sans", linespacing=1.55)

    # ---- sliders (bottom-left) --------------------------------------------------------
    def _build_sliders(self):
        from matplotlib.widgets import Slider
        self.sliders = {}
        y = 0.37
        for (key, label, vmin, vmax, step, level, _m) in PARAMS:
            sax = self.fig.add_axes([0.14, y, 0.33, 0.02])
            s = Slider(sax, label, vmin, vmax, valinit=self.vals[key], valstep=step)
            s.label.set_fontsize(9)
            s.valtext.set_fontsize(9)
            s.on_changed(lambda val, k=key, lv=level: self._on_slider(k, lv, val))
            self.sliders[key] = s
            y -= 0.037

    # ---- buttons (bottom-right, under the meanings) ----------------------------------
    def _build_buttons(self):
        from matplotlib.widgets import Button
        self.buttons = {}
        specs = [
            ("save", "Save to config", [0.55, 0.115, 0.21, 0.055], self._save_config),
            ("reset", "Reset params", [0.78, 0.115, 0.21, 0.055], self._reset),
            ("prev", "◀ Prev event", [0.55, 0.05, 0.21, 0.055], lambda e: self._nav(-1)),
            ("next", "Next event ▶", [0.78, 0.05, 0.21, 0.055], lambda e: self._nav(+1)),
        ]
        for name, txt, rect, cb in specs:
            b = Button(self.fig.add_axes(rect), txt)
            b.label.set_fontsize(10)
            b.on_clicked(cb)
            self.buttons[name] = b

    # ---- data + compute ---------------------------------------------------------------
    def _load_window(self):
        cfg = self.cfg
        ev = self.cands[self.idx]
        t_disp0 = ev["start_time"] - timedelta(seconds=_PAD_PRE)
        t_disp1 = ev["end_time"] + timedelta(seconds=_PAD_POST)
        t_load0 = t_disp0 - timedelta(seconds=_LEAD)
        Dfull, t_actual0 = read_window(t_load0, t_disp1, cfg)
        self.Dfull = Dfull
        self.n_lead = int(round((t_disp0 - t_actual0).total_seconds() * self.fs))
        self.n_lead = min(max(self.n_lead, 0), max(0, Dfull.shape[0] - 1))
        self.nch = Dfull.shape[1]
        # raw-waveform + envelope backgrounds (independent of the STA/LTA params) so you can look
        # at the actual signal — not only the CF — while judging how a parameter change moves things.
        self.D = Dfull[self.n_lead:, :]
        sc = np.percentile(np.abs(self.D), 95, axis=0) + 1e-9
        self.amp = self.D / sc
        env = np.abs(hilbert(self.D, axis=0))
        self.env = env / (env.max(axis=0, keepdims=True) + 1e-9)

    def _recompute(self, level):
        v = self.vals
        if level == "cf":
            nsta = max(1, int(v["sta_len"] * self.fs))
            nlta = max(2, int(v["lta_len"] * self.fs))
            csl_full = sta_lta_centred_multichan(self.Dfull, nsta, nlta)
            self.csl = csl_full[self.n_lead:, :]
            level = "track"
        if level == "track":
            tp = {"neighbor_tol": int(v["neighbor_tol"]),
                  "time_tol": float(v["time_tol"]),
                  "slope_tol": float(v["slope_tol"])}
            mp = {"merge_time_tol": float(v["merge_time_tol"])}
            out = trig_allchan_and_collapsed_index_and_save_wavefonts(
                self.csl, float(v["threshold_level"]), float(v["minimal_event_interval"]),
                self.fs, self.dx, tp, mp)
            self.wavefronts = out[-1]
            self.empty_chan = out[-2]
        self._draw()

    # ---- drawing ----------------------------------------------------------------------
    def _draw(self):
        ax = self.ax
        ax.clear()
        nt = self.csl.shape[0]
        t1 = nt / self.fs
        dist = (self.d1 + np.arange(self.nch) * self.dx) / 1000.0

        # background: 1 raw amplitude · 2 STA/LTA CF · 3 envelope · 4 wiggle traces.
        # `-`/`=` scale the intensity: image bg -> tighter/looser colour range; wiggle -> trace gain.
        g = self.gain[self.bg[0]]
        if self.bg[0] == "4":
            nm = f"wiggle ({min(self.n_wig, self.nch)} traces · gain ×{g:.2f})"
            self._draw_wiggle_bg(ax, t1, dist, g)
        else:
            if self.bg[0] == "1":
                img, cmap, nm = self.amp.T, "RdBu_r", "amplitude (raw)"
                vmn, vmx = -0.6 / g, 0.6 / g
            elif self.bg[0] == "3":
                img, cmap, nm = self.env.T, "magma", "envelope"
                vmn, vmx = 0, 1.0 / g
            else:
                cf = np.clip(self.csl / (np.percentile(self.csl, 99, axis=0, keepdims=True) + 1e-9), 0, 1)
                img, cmap, nm = cf.T, "magma", "STA/LTA CF"
                vmn, vmx = 0, 1.0 / g
            nm = f"{nm} (contrast ×{g:.2f})"
            ax.imshow(img, aspect="auto", cmap=cmap, vmin=vmn, vmax=vmx,
                      extent=[0, t1, dist.max(), dist.min()], interpolation="nearest")

        # Overlay the tracked wavefronts FINELY and SEMI-TRANSPARENTLY, so the per-channel picks and
        # the wavefront line never mask the (very detailed) waveform underneath. Each wavefront = a
        # thin line through its per-channel picks + small dots. lime = passes min_channels, grey = below.
        mc = self.vals["min_channels"]
        npass = 0
        for wf in self.wavefronts:
            det = np.asarray(wf["details"], dtype=float)
            if det.size == 0:
                continue
            det = det[np.argsort(det[:, 0])]
            dk = (self.d1 + det[:, 0] * self.dx) / 1000.0
            tt = det[:, 1] / self.fs
            ok = wf["score"] >= mc
            npass += int(ok)
            col = "lime" if ok else "0.75"
            ax.plot(tt, dk, color=col, lw=1.1 if ok else 0.6,
                    alpha=0.55 if ok else 0.3, zorder=4 if ok else 3)
            ax.scatter(tt, dk, s=5 if ok else 3, c=col, alpha=0.5 if ok else 0.28,
                       edgecolors="none", zorder=4 if ok else 3)

        ev = self.cands[self.idx]
        ax.set_xlim(0, t1)
        ax.set_ylim(dist.max(), dist.min())
        ax.set_xlabel("time in window (s)")
        ax.set_ylabel("distance along cable (km)")
        self.status = (f"tracked {len(self.wavefronts)} wavefronts · "
                       f"{npass} pass min_channels(≥{int(mc)}) · "
                       f"empty {getattr(self, 'empty_chan', 0)}/{self.nch} ch")
        ax.set_title(
            f"TUNE  ·  event {self.idx+1}/{len(self.cands)}  "
            f"{ev['start_time'].strftime('%Y-%m-%d %H:%M:%S')}  ·  bg: {nm}   ·   press  h  for keys\n"
            f"{self.status}   (lime = passes min_channels · grey = tracked but below it)",
            fontsize=11)
        if self.help_visible[0]:
            draw_shortcut_overlay(ax, _TUNE_HELP_ROWS)
        self.fig.canvas.draw_idle()

    def _draw_wiggle_bg(self, ax, t1, dist, gmul=1.0):
        """Wiggle background: sampled raw traces at their true distance (density = , / . ; amplitude
        gain = - / =), so you can see the individual waveforms the detector triggers on. The base
        gain scales to the inter-trace spacing so traces don't overlap at ×1."""
        D = self.D
        n = min(self.n_wig, self.nch)
        ksel = np.unique(np.linspace(0, self.nch - 1, max(n, 1)).astype(int))
        dk = (self.d1 + ksel * self.dx) / 1000.0
        t = np.arange(D.shape[0]) / self.fs
        base = 0.45 * (dk.max() - dk.min()) / max(len(dk) - 1, 1) if len(dk) > 1 else 0.45
        gain = base * gmul
        for ki, d0 in zip(ksel, dk):
            tr = D[:, ki].astype(float)
            tr = tr / (np.max(np.abs(tr)) + 1e-9) * gain
            ax.plot(t, d0 + tr, "k", lw=0.4, zorder=0.5)
            ax.fill_between(t, d0, d0 + tr, where=(tr > 0), color="k", linewidth=0, zorder=0.5)

    # ---- callbacks --------------------------------------------------------------------
    def _on_slider(self, key, level, val):
        self.vals[key] = float(val)
        if level == "filter":
            self._draw()                 # min_channels: just re-colour, no recompute
        else:
            self._recompute(level)

    def _nav(self, step):
        self.idx = int(np.clip(self.idx + step, 0, len(self.cands) - 1))
        self._load_window()
        self._recompute("cf")

    def _reset(self, _event=None):
        for (key, *_r) in PARAMS:
            self.vals[key] = float(_cfg_value(self.cfg, key))
            if key in self.sliders:
                self.sliders[key].set_val(self.vals[key])   # fires _on_slider -> recompute
        self._recompute("cf")

    def _save_config(self, _event=None):
        path = self.cfg.config_path
        data = {}
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
        v = self.vals
        data["threshold_level"] = float(v["threshold_level"])
        data["sta_len"] = float(v["sta_len"])
        data["lta_len"] = float(v["lta_len"])
        data["minimal_event_interval"] = float(v["minimal_event_interval"])
        data["min_channels"] = int(v["min_channels"])
        data["tracking_params"] = {"neighbor_tol": int(v["neighbor_tol"]),
                                   "time_tol": float(v["time_tol"]),
                                   "slope_tol": float(v["slope_tol"])}
        data["merge_params"] = {"merge_time_tol": float(v["merge_time_tol"])}
        with open(path, "w") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        msg = f"saved parameters to {path}"
        print(msg)
        self.ax.set_title(self.ax.get_title().split("\n")[0] + f"\n✓ {msg}", fontsize=11)
        self.fig.canvas.draw_idle()

    def run(self):
        print(__doc__)
        print(f">>> tuning on {len(self.cands)} candidate windows · Save writes -> {self.cfg.config_path}")
        self.plt.show()


def selftest(cfg, cands, out_png, start_idx=0):
    """Headless: build the tuner, exercise a couple of parameter changes, render a PNG."""
    t = Tuner(cfg, cands, start_idx=start_idx, interactive=False)
    n0 = len(t.wavefronts)
    t._on_slider("threshold_level", "track", t.vals["threshold_level"] + 20)
    n1 = len(t.wavefronts)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    t.fig.savefig(out_png, dpi=110)
    return {"png": out_png, "nch": t.nch, "wavefronts_default": n0,
            "wavefronts_higher_threshold": n1, "status": t.status}
