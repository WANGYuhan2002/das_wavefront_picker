"""
Interactive wavefront picker — ONE window, ONE wavefront, ONE label, ONE file.

The atomic unit here is the WAVEFRONT, not an "event". Detection produces candidate
wavefronts; for each you open its raw-DAS window, refine the single wavefront by clicking /
dragging anchors, assign a LABEL (P / S / vehicle / explosive / noise / unknown — from a
controlled vocabulary), and save it as its own `.h5` file. There is deliberately no P/S
pairing and no event grouping — associating wavefronts is a separate, downstream problem.

The window seeds its wavefront from the Method-B tracked wavefront if present, else from a
per-channel STA/LTA auto-pick. Four switchable backgrounds share the distance-time axes:
1 phase-adaptive amplitude, 2 STA-LTA CF, 3 envelope (all image, with an extra `lta_len`
lead-in so the CF is warmed up), and 4 a wiggle gather (sampled traces at their true distance,
with the wavefront overlaid on top). Local (partial-array) wavefronts auto-focus onto their span.

Runs interactively (native matplotlib window) or headless (`selftest`) for verification.
"""

import os
import copy
import hashlib
import numpy as np
from datetime import timedelta

import matplotlib
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import median_filter

from stalta import sta_lta_centred          # vendored
from . import catalog as cat
from .dasio import read_window


# ============================================================= raw-DAS window =========
class WaveWindow:
    """Raw-DAS window for one candidate wavefront + its backgrounds + a seed wavefront."""

    def __init__(self, cand, cfg, exp_dir):
        self.cand = cand
        self.cfg = cfg
        fs = cfg.fs
        lead_s = cfg.lta_len                          # CF warm-up lead

        t_disp0 = cand["start_time"] - timedelta(seconds=cfg.pad_pre_s)
        t_disp1 = cand["end_time"] + timedelta(seconds=cfg.pad_post_s)
        t_load0 = t_disp0 - timedelta(seconds=lead_s)

        # read with the CF warm-up lead; align to the TRUE returned start time so a wavefront
        # near a file/gap boundary is not silently time-shifted (the loader does not front-pad).
        Dfull, t_actual0 = read_window(t_load0, t_disp1, cfg)
        n_lead = int(round((t_disp0 - t_actual0).total_seconds() * fs))
        n_lead = min(max(n_lead, 0), max(0, Dfull.shape[0] - 1))

        self.D = Dfull[n_lead:, :]                     # display window [nt, nch]
        self.nt, self.nch = self.D.shape
        self.fs = fs
        disp_start = t_actual0 + timedelta(seconds=n_lead / fs)   # real display start
        self.t = np.arange(self.nt) / fs               # display time (s), 0 = disp_start
        # candidate reference time sits at t = ref_mark (== cfg.pad_pre_s except at a boundary)
        self.ref_mark = (cand["start_time"] - disp_start).total_seconds()
        self.t_disp0 = disp_start
        self.chan = np.arange(self.nch)
        self.dist_km = (cfg.d1 + self.chan * cfg.dx) / 1000.0

        # --- characteristic function (raw MAD-normalised) on the full load, then trim ---
        nsta, nlta = cfg.nsta, cfg.nlta
        CFfull = np.empty_like(Dfull, dtype=float)
        for k in range(self.nch):
            CFfull[:, k] = sta_lta_centred(Dfull[:, k], nsta, nlta)
        self.CF_raw = CFfull[n_lead:, :]               # raw CF, [nt, nch]
        self.CF_disp = self.CF_raw / (np.percentile(self.CF_raw, 99, axis=0, keepdims=True) + 1e-9)
        self.CF_disp = np.clip(self.CF_disp, 0, 1)
        env = np.abs(_hilbert_mag(self.D))
        self.ENV = env / (env.max(axis=0, keepdims=True) + 1e-9)

        # --- single seed wavefront (display seconds per channel) --------------------
        b_seed = cat.read_seed_wavefront(exp_dir, cand)
        rel = cat.seed_to_relative(b_seed, fs) if b_seed is not None else None
        if rel is not None:
            self.seed, self.seed_conf = self._seed_from_tracked(rel)
        else:
            self.seed, self.seed_conf = self._seed_auto()

    # -- seeds --------------------------------------------------------------------------
    def _seed_from_tracked(self, rel):
        """rel = (channels, arrival_rel_to_ref_s, values). Densify across all channels."""
        ch, arr_rel, _val = rel
        arr_disp = self.ref_mark + arr_rel                     # into display coords
        ch = np.clip(ch, 0, self.nch - 1)
        # flat-extrapolated interpolation over channel index (safe: clipped to endpoints)
        dense = np.interp(self.chan, ch, arr_disp, left=arr_disp[0], right=arr_disp[-1])
        conf = (self.chan >= ch.min()) & (self.chan <= ch.max())
        return np.clip(dense, self.t[0], self.t[-1]), conf

    def _seed_auto(self):
        """Per-channel STA/LTA argmax pick within the display window."""
        arr = np.empty(self.nch)
        val = np.empty(self.nch)
        for k in range(self.nch):
            j = int(np.argmax(self.CF_raw[:, k]))
            arr[k] = self.t[j]
            val[k] = self.CF_raw[j, k]
        arr = median_filter(arr, 7)
        conf = val > 0.5 * self.cfg.threshold_level
        return np.clip(arr, self.t[0], self.t[-1]), conf

    # -- amplitude background (arrival-adaptive) ----------------------------------------
    def amp_background(self, centre_t=None):
        """Per-channel amplitude normalised to the 95th pct of |amp| in a region around
        `centre_t` (the wavefront's mean pick), so the arrival saturates into colour."""
        if centre_t is None:
            m = np.ones(self.nt, bool)
        else:
            m = (self.t >= centre_t - 0.6) & (self.t <= centre_t + 1.6)
            if m.sum() < 5:
                m = np.ones(self.nt, bool)
        sc = np.percentile(np.abs(self.D[m, :]), 95, axis=0) + 1e-9
        return self.D / sc


def _hilbert_mag(x):
    from scipy.signal import hilbert
    return hilbert(x, axis=0)


def _resolve_out_dir(cfg, exp_dir):
    """Where to WRITE picks. Beside the catalog if it is inside our own Output tree;
    otherwise a unique, writable folder under the app Output keyed to the source path (so we
    never write into an external / read-only catalog folder and never collide across sources)."""
    out_root = os.path.abspath(cfg.output_dir)
    exp_abs = os.path.abspath(exp_dir)
    if exp_abs == out_root or exp_abs.startswith(out_root + os.sep):
        return exp_dir
    tag = os.path.basename(os.path.normpath(exp_dir)) or "external"
    h = hashlib.md5(exp_abs.encode()).hexdigest()[:6]
    return os.path.join(cfg.output_dir, f"picks_{tag}_{h}")


_HELP_ROWS = [
    ("Navigate",   "n / b",     "next / previous wavefront   (auto-saves if edited)"),
    ("Set pick",   "click",     "pull the wavefront to the click (nearest anchor)"),
    ("Anchors",    "a / x",     "add / delete anchor at the mouse"),
    ("Undo/Reset", "u / r",     "undo   ·   reset to the auto seed"),
    ("Label",      "l",         "cycle label  (P · S · vehicle · explosive · noise · unknown)"),
    ("Background", "1 2 3 4",   "amplitude · STA-LTA CF · envelope · wiggle traces"),
    ("Intensity",  "- / =",     "dim / brighten active bg  (colour range · wiggle gain)"),
    ("Wiggle #",   ", / .",     "wiggle trace density: fewer / more traces"),
    ("View",       "d / f",     "toggle seed line  ·  focus (local wavefront)"),
    ("Confidence", "c",         "toggle low-confidence  ±0.15 km  at the mouse"),
    ("File",       "w / q",     "SAVE this wavefront   ·   quit"),
]

# Fold full-width punctuation/digits (produced when a CJK IME is on) to their ASCII keys, so the
# shortcuts still work without switching the input method to English.
_FULLWIDTH = {"，": ",", "。": ".", "．": ".", "、": ".", "－": "-", "—": "-", "＝": "=",
              "＋": "+", "１": "1", "２": "2", "３": "3", "４": "4"}


def draw_shortcut_overlay(ax, rows, title="KEYBOARD SHORTCUTS"):
    """Grouped shortcut overlay drawn on `ax` (axes-fraction coords). Columns are positioned by
    x (not space-padding), so alignment is font-independent. `rows` = [(label, keys, desc), ...].
    Shared by the picker (`h`) and the tuner (`h`)."""
    from matplotlib.patches import FancyBboxPatch
    x0, x1, y0, y1 = 0.15, 0.85, 0.17, 0.83
    ax.add_patch(FancyBboxPatch(
        (x0, y0), x1 - x0, y1 - y0, transform=ax.transAxes, zorder=20,
        boxstyle="round,pad=0.008,rounding_size=0.02",
        facecolor="white", alpha=0.95, edgecolor="0.35", linewidth=1.2))
    xl, xk, xd = x0 + 0.03, x0 + 0.17, x0 + 0.29
    ax.text(xl, y1 - 0.045, title, transform=ax.transAxes,
            ha="left", va="center", fontsize=12, fontweight="bold", zorder=21)
    ax.text(x1 - 0.03, y1 - 0.045, "press  h  to close", transform=ax.transAxes,
            ha="right", va="center", fontsize=9.5, color="0.45", zorder=21)
    ytop = y1 - 0.115
    kw = dict(transform=ax.transAxes, va="top", ha="left", ma="left",
              fontsize=10.5, zorder=21, linespacing=1.75)
    ax.text(xl, ytop, "\n".join(r[0] for r in rows), fontweight="bold", **kw)
    ax.text(xk, ytop, "\n".join(r[1] for r in rows), family="monospace", **kw)
    ax.text(xd, ytop, "\n".join(r[2] for r in rows), **kw)

_WF_COLOR = "lime"


# ============================================================= the picker =============
class Picker:
    def __init__(self, cfg, cands, exp_dir, start_idx=0, interactive=True):
        self.cfg = cfg
        self.cands = cands                            # list of candidate wavefronts to browse
        self.exp_dir = exp_dir                        # SEEDs are read from here (read-only OK)
        # Picks are written beside the catalog when it is inside our own Output tree; for an
        # EXTERNAL catalog (e.g. a read-only matched_filter_release folder) we never write
        # there — output goes to a unique, writable app-Output folder keyed to the source path.
        self.out_dir = _resolve_out_dir(cfg, exp_dir)
        self.picks_dir = os.path.join(self.out_dir, "Picks")
        self.index_path = os.path.join(self.out_dir, "picks_index.csv")
        self.interactive = interactive
        self.labels = list(cfg.labels)
        self.label = [cfg.default_label if cfg.default_label in self.labels else self.labels[0]]
        self.bg = ["1"]
        self.mouse = [None]
        self.show_seed = [True]
        self.help_visible = [False]
        self._cache = {}
        self.idx = None
        self.win = None
        self.state = {}
        self.hist = []
        self.conf = None
        self.dirty = False

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
        # main axes + a slim colorbar axes (fixed positions so nothing shifts on redraw)
        self.fig = plt.figure(figsize=(15, 9))
        self.ax = self.fig.add_axes([0.06, 0.09, 0.86, 0.80])
        self.cax = self.fig.add_axes([0.935, 0.09, 0.012, 0.80])
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        self._sm = ScalarMappable(norm=Normalize(0, 1), cmap="magma")
        self._sm.set_array([])
        self.cbar = self.fig.colorbar(self._sm, cax=self.cax)
        self.cbar.ax.tick_params(labelsize=8)
        # display tunables (persist across wavefronts): per-background intensity/gain + wiggle count
        self.gain = {"1": 1.0, "2": 1.0, "3": 1.0, "4": 1.0}
        self.n_wig = int(cfg.n_wiggle)
        self._n_shown = 0                              # wiggle traces visible in the current view

        self.load(start_idx)
        if interactive:
            self.fig.canvas.mpl_connect("motion_notify_event", self.on_move)
            self.fig.canvas.mpl_connect("button_press_event", self.on_click)
            self.fig.canvas.mpl_connect("key_press_event", self.on_key)

    # -- window loading -----------------------------------------------------------------
    def load(self, idx):
        idx = int(np.clip(idx, 0, len(self.cands) - 1))
        if idx not in self._cache:
            if len(self._cache) > 4:                 # keep the cache small
                self._cache.pop(next(iter(self._cache)))
            self._cache[idx] = WaveWindow(self.cands[idx], self.cfg, self.exp_dir)
        self.idx = idx
        self.win = self._cache[idx]
        dist = self.win.dist_km

        # LOCAL wavefront handling: if it only lights up part of the array, concentrate the
        # anchors on the supported span (+ a small margin) instead of spreading them uniformly
        # over the whole array, and auto-focus the view onto that span.
        conf0 = self.win.seed_conf
        na = self.cfg.n_anchors
        self.support = None
        self.is_local = False
        if conf0.any() and conf0.sum() < 0.9 * len(conf0):
            self.is_local = True
            lo, hi = dist[conf0].min(), dist[conf0].max()
            step = dist[1] - dist[0] if len(dist) > 1 else 0.0
            margin = max(0.08 * (hi - lo), 3 * step)
            lo2, hi2 = max(dist.min(), lo - margin), min(dist.max(), hi + margin)
            self.support = (lo2, hi2)
            inner = np.linspace(lo2, hi2, max(na - 2, 2))
            anch_s0 = np.unique(np.concatenate([[dist.min()], inner, [dist.max()]]))
            self.focus = [True]
        else:
            anch_s0 = np.linspace(dist.min(), dist.max(), na)
            self.focus = [False]

        self.state = {"s": anch_s0.copy(), "t": np.interp(anch_s0, dist, self.win.seed)}
        self.auto = copy.deepcopy(self.state)
        self.hist = []
        self.conf = conf0.copy()
        self.dirty = False
        self.redraw()

    # -- wavefront math -----------------------------------------------------------------
    def wavefront_of(self, a_s, a_t):
        u, ui = np.unique(a_s, return_index=True)
        if len(u) < 2:
            return np.full_like(self.win.dist_km, a_t[0])
        f = PchipInterpolator(u, a_t[ui], extrapolate=True)
        wf = median_filter(f(self.win.dist_km), 7)
        return np.clip(wf, self.win.t[0], self.win.t[-1])

    def wavefront(self):
        return self.wavefront_of(self.state["s"], self.state["t"])

    def push(self):
        self.hist.append(copy.deepcopy(self.state))
        self.dirty = True

    # -- drawing ------------------------------------------------------------------------
    def redraw(self):
        ax = self.ax
        w = self.win
        ax.clear()
        g = self.gain[self.bg[0]]                      # per-background intensity/gain
        if self.bg[0] == "4":                          # wiggle gather (sampled traces)
            nm = "wiggle (sampled traces)"
            self._draw_wiggle_bg(ax, g)
            self.cax.set_visible(False)                # no colour scale for the wiggle view
        else:                                          # image backgrounds 1 / 2 / 3
            if self.bg[0] == "1":
                centre = float(np.mean(self.wavefront()))
                img, cmap, base, nm = w.amp_background(centre).T, "RdBu_r", 0.6, \
                    "amplitude (arrival-scaled)"
                vmn, vmx = -base / g, base / g         # higher gain -> tighter range -> more contrast
            elif self.bg[0] == "2":
                img, cmap, nm = w.CF_disp.T, "magma", "STA/LTA CF"
                vmn, vmx = 0, 1.0 / g
            else:
                img, cmap, nm = w.ENV.T, "magma", "envelope"
                vmn, vmx = 0, 1.0 / g
            ax.imshow(img, aspect="auto", cmap=cmap, vmin=vmn, vmax=vmx,
                      extent=[w.t[0], w.t[-1], w.dist_km.max(), w.dist_km.min()],
                      interpolation="nearest")
            self._sm.set_cmap(cmap)                     # keep the colorbar in sync with the image
            self._sm.set_clim(vmn, vmx)
            self.cbar.update_normal(self._sm)
            self.cax.set_visible(True)

        # low-confidence distance bands (grey)
        cs = self.conf
        edges = np.where(np.diff(cs.astype(int)) != 0)[0]
        bnds = np.concatenate([[0], edges + 1, [len(cs)]])
        for i in range(len(bnds) - 1):
            lo, hi = bnds[i], bnds[i + 1] - 1
            if not cs[lo]:
                ax.axhspan(w.dist_km[lo], w.dist_km[hi], color="0.5", alpha=0.18, zorder=1)

        ax.axvline(w.ref_mark, color="w", ls=":", lw=0.8, alpha=0.6)   # candidate ref time

        wa = float(getattr(self.cfg, "wavefront_alpha", 0.6))   # line opacity (config-tunable)
        wf = self.wavefront()
        # dashed reference = the wavefront through the INITIAL anchors, i.e. exactly what `r`
        # resets to (NOT the raw per-channel seed, which the 16-anchor model can't reproduce) —
        # so after reset the solid line lies right on the dashed line.
        if self.show_seed[0]:
            seed_wf = self.wavefront_of(self.auto["s"], self.auto["t"])
            ax.plot(seed_wf, w.dist_km, color=_WF_COLOR, ls="--", lw=0.9, alpha=min(wa, 0.5),
                    zorder=3)
        ax.plot(wf, w.dist_km, color=_WF_COLOR, lw=2.5, alpha=wa, zorder=4,
                label=f"wavefront · {self.label[0]}")
        ax.scatter(self.state["t"], self.state["s"], s=45, c=_WF_COLOR, edgecolor="k",
                   alpha=min(1.0, wa + 0.3), zorder=5)

        cand = self.cands[self.idx]
        ax.set_xlim(w.t[0], w.t[-1])
        if self.focus[0] and self.support is not None:
            lo, hi = self.support
            ax.set_ylim(hi, lo)                        # zoom to the local wavefront's span
        else:
            ax.set_ylim(w.dist_km.max(), w.dist_km.min())
        ax.set_xlabel("time in window (s)   ·   dotted = candidate reference time")
        ax.set_ylabel("distance along cable (km)")
        tstr = cand["start_time"].strftime("%Y-%m-%d %H:%M:%S")

        # Title carries STATE only (identity + label + view); shortcuts live in the h overlay.
        parts = [f"WAVEFRONT {self.idx+1}/{len(self.cands)}", tstr,
                 f"{cand['source']}·{cand['n_channels']}ch"]
        if self.is_local:
            parts.append(f"LOCAL[{'focus' if self.focus[0] else 'full'}]")
        if self.dirty:
            parts.append("UNSAVED")
        line1 = "   ·   ".join(parts)
        if self.bg[0] == "4":
            tune = f"gain ×{g:.2f} · {self._n_shown} traces shown"
        else:
            tune = f"contrast ×{g:.2f}"
        line2 = (f"label: {self.label[0]}   ·   bg: {nm} ({tune})"
                 "          press  h  for shortcuts")
        ax.set_title(line1 + "\n" + line2, fontsize=11,
                     color="firebrick" if self.dirty else "black")
        ax.legend(loc="upper right", fontsize=10)

        if self.help_visible[0]:
            self._draw_help(ax)
        self.fig.canvas.draw_idle()

    def _draw_help(self, ax):
        draw_shortcut_overlay(ax, _HELP_ROWS)

    def _wiggle_krange(self):
        """Channel-index range [k0, k1] the wiggle background samples over: ALWAYS the whole
        array, so `,`/`.` are a plain density knob independent of the local-focus state (focus
        only zooms the y-view). (k1 - k0 + 1) is the max # of traces available."""
        return 0, self.win.nch - 1

    def _draw_wiggle_bg(self, ax, gmul=1.0):
        """Background '4': sampled DAS traces drawn at their TRUE distance (wiggle gather), in
        the same distance-time axes as the image backgrounds, so the wavefront overlays on top
        and crosses each trace at the arrival. `,`/`.` set the trace density (`self.n_wig`) over
        the whole array and `gmul` (via `-`/`=`) the amplitude; the base gain scales to the
        inter-trace distance so traces don't overlap at ×1. Focus just zooms the visible range."""
        w = self.win
        k0, k1 = self._wiggle_krange()
        n = min(self.n_wig, k1 - k0 + 1)
        ksel = np.unique(np.linspace(k0, k1, max(n, 1)).astype(int))
        d = w.dist_km[ksel]
        # amplitude scales to the inter-trace spacing WITHIN THE VISIBLE range (not the whole
        # array), so traces are well-sized whether you're zoomed into a local wavefront or not.
        if self.focus[0] and self.support is not None:
            vlo, vhi = self.support
        else:
            vlo, vhi = w.dist_km.min(), w.dist_km.max()
        n_show = int(np.sum((d >= vlo) & (d <= vhi)))
        self._n_shown = n_show                          # traces actually visible (for the title)
        gain = 0.45 * (vhi - vlo) / max(n_show, 1) * gmul
        for ki, dk in zip(ksel, d):
            tr = w.D[:, ki].astype(float)
            tr = tr / (np.max(np.abs(tr)) + 1e-9) * gain
            ax.plot(w.t, dk + tr, "k", lw=0.4, zorder=0.5)
            ax.fill_between(w.t, dk, dk + tr, where=(tr > 0), color="k", linewidth=0, zorder=0.5)

    # -- interaction --------------------------------------------------------------------
    def on_move(self, ev):
        if ev.inaxes == self.ax and ev.ydata is not None:
            self.mouse[0] = ev.ydata

    def on_click(self, ev):
        if ev.inaxes != self.ax or ev.xdata is None or ev.button != 1:
            return
        a = self.state
        k = int(np.argmin(np.abs(a["s"] - ev.ydata)))
        self.push()
        a["t"][k] = ev.xdata
        self.redraw()

    def on_key(self, ev):
        a = self.state
        ms = ev.ydata if (ev.inaxes == self.ax and ev.ydata is not None) else self.mouse[0]
        k = _FULLWIDTH.get(ev.key, ev.key)             # tolerate a CJK IME's full-width keys
        if k in ("1", "2", "3", "4"):
            self.bg[0] = k
            self.redraw()
        elif k in ("=", "+", "-"):                     # intensity/gain of the active background
            f = 1.3 if k in ("=", "+") else 1 / 1.3
            self.gain[self.bg[0]] = float(np.clip(self.gain[self.bg[0]] * f, 0.1, 12.0))
            self.redraw()
        elif k in (",", "."):                          # wiggle-background trace count
            k0, k1 = self._wiggle_krange()
            avail = k1 - k0 + 1                         # can't show more traces than channels in view
            f = 1.4 if k == "." else 1 / 1.4
            self.n_wig = int(np.clip(round(self.n_wig * f), 4, avail))
            if self.bg[0] == "4":
                self.redraw()
        elif k == "l":
            i = self.labels.index(self.label[0]) if self.label[0] in self.labels else -1
            self.label[0] = self.labels[(i + 1) % len(self.labels)]
            self.dirty = True
            self.redraw()
        elif k == "d":
            self.show_seed[0] = not self.show_seed[0]
            self.redraw()
        elif k == "f":
            if self.support is not None:
                self.focus[0] = not self.focus[0]
                self.redraw()
        elif k == "r":
            self.push()
            self.state = copy.deepcopy(self.auto)
            self.redraw()
        elif k == "u":
            if self.hist:
                self.state = self.hist.pop()
                self.redraw()
        elif k == "a":
            if ms is None:
                return
            nt = float(np.interp(ms, self.win.dist_km, self.wavefront()))
            self.push()
            a["s"] = np.append(a["s"], ms)
            a["t"] = np.append(a["t"], nt)
            self.redraw()
        elif k == "x":
            if ms is None or len(a["s"]) <= 2:
                return
            j = int(np.argmin(np.abs(a["s"] - ms)))
            self.push()
            a["s"] = np.delete(a["s"], j)
            a["t"] = np.delete(a["t"], j)
            self.redraw()
        elif k == "c":
            if ms is None:
                return
            band = np.abs(self.win.dist_km - ms) <= 0.15
            self.conf[band] = ~self.conf[band].all()
            self.dirty = True
            self.redraw()
        elif k == "w":
            self.save()
        elif k in ("n", "b", "[", "]"):
            if self.dirty:
                self.save()
            step = 1 if k in ("n", "]") else -1
            self.load(self.idx + step)
        elif k == "h":
            self.help_visible[0] = not self.help_visible[0]
            self.redraw()
        elif k == "q":
            self.plt.close(self.fig)

    # -- export -------------------------------------------------------------------------
    def save(self):
        """Persist THIS wavefront as its own `.h5` file (release-compatible `wavefront_data`
        + arrival_rel_s + conf + metadata incl. the label) and update the picks index CSV."""
        cfg = self.cfg
        cand = self.cands[self.idx]
        w = self.win
        ref_time = cand["start_time"]
        label = self.label[0]

        wf = self.wavefront()                              # display seconds per channel
        arr_rel_s = wf - w.ref_mark                         # relative to candidate ref time
        samp = np.clip((wf * cfg.fs).astype(int), 0, w.nt - 1)
        vals = w.CF_raw[samp, np.arange(w.nch)]             # CF value at each pick
        keep = self.conf
        n_conf = int(keep.sum())

        # apparent slowness along the cable (s/km) + mean arrival, over confident channels
        if n_conf >= 2:
            slope = float(np.polyfit(w.dist_km[keep], arr_rel_s[keep], 1)[0])
            mean_arr = float(np.mean(arr_rel_s[keep]))
        else:
            slope = float("nan")
            mean_arr = float(np.mean(arr_rel_s))

        meta = {"source": cand["source"], "app_slowness_s_per_km": slope,
                "mean_arrival_s": mean_arr, "n_conf": n_conf}
        fname = cat.write_wavefront_h5(self.picks_dir, ref_time, label, w.chan,
                                       arr_rel_s, vals, cfg.fs, conf=self.conf, meta=meta)
        cat.update_picks_index(self.index_path, {
            "file": fname, "ref_time": ref_time.isoformat(), "label": label,
            "source": cand["source"], "n_channels": int(w.nch), "n_conf": n_conf,
            "mean_arrival_s": mean_arr, "app_slowness_s_per_km": slope,
        })
        fig_path = os.path.join(self.picks_dir, os.path.splitext(fname)[0] + ".png")
        self.fig.savefig(fig_path, dpi=120)
        self.dirty = False
        print(f"SAVED wavefront {self.idx+1} -> {fname}  "
              f"label={label}  ({n_conf}/{w.nch} ch reliable)")
        if self.interactive:
            self.redraw()
        return fname

    def run(self):
        print(__doc__)
        print(f">>> {len(self.cands)} candidate wavefronts · picks -> {self.picks_dir}")
        self.plt.show()


# ============================================================= headless selftest ======
def selftest(cfg, cands, exp_dir, out_png, idx=0):
    """Non-interactive smoke test: load one candidate, build backgrounds + seed, set a label,
    render a PNG, and export the wavefront file. Exercises the whole pick/export path."""
    p = Picker(cfg, cands, exp_dir, start_idx=idx, interactive=False)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    p.fig.savefig(out_png, dpi=120)
    fname = p.save()
    return {"png": out_png, "file": fname, "label": p.label[0],
            "nt": p.win.nt, "nch": p.win.nch,
            "cf_max": float(np.nanmax(p.win.CF_raw))}
