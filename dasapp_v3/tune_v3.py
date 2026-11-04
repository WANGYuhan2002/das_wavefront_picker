"""
v3 parameter tuner = the v2 Tuner + a live channel-quality mask toggle.

Additions (dasapp/tune.py is NOT modified):
  * the resolved channel mask is applied to the CF before tracking, exactly as in
    v3 detection (csl[:, masked] = 0), so the preview shows what `run_v3.py detect`
    will actually do;
  * press `m` to toggle the mask ON/OFF and watch the tracking result change — the
    single most direct way to see what the mask buys (or costs) on a given event;
  * masked channels show as red stripes at the left edge of the waterfall; the status
    line reports the mask state.

Two v2 methods are copied with marked edits (`# V3`): _recompute (mask injection —
the v2 method computes and consumes the CF in one flow with no hook) and _draw
(stripes + status + help rows). Everything else is inherited untouched.
"""

import numpy as np

from stalta import (                                    # vendored (v2 path)
    sta_lta_centred_multichan,
    trig_allchan_and_collapsed_index_and_save_wavefonts,
)
from dasapp.tune import Tuner, _TUNE_HELP_ROWS
from dasapp.picker import draw_shortcut_overlay

_V3_HELP_ROWS = list(_TUNE_HELP_ROWS) + [
    ("Mask", "m", "toggle the channel-quality mask on/off (v3)"),
]


class TunerV3(Tuner):

    def __init__(self, cfg, cands, start_idx=0, interactive=True):
        mask = getattr(cfg, "chan_mask", None)
        self.mask_full = None if mask is None else np.asarray(mask, bool)
        self.mask_on = self.mask_full is not None and self.mask_full.any()
        self.csl_raw = None
        super().__init__(cfg, cands, start_idx=start_idx, interactive=interactive)

    # -- mask helpers -------------------------------------------------------------------
    def _mask_vec(self):
        """Mask sliced to this window's channel count (or None)."""
        if self.mask_full is None:
            return None
        return self.mask_full[: self.nch]

    # V3-copy of Tuner._recompute (dasapp/tune.py) — identical flow, plus the CF is kept
    # unmasked in self.csl_raw and the mask is applied into self.csl when mask_on.
    def _recompute(self, level):
        v = self.vals
        if level == "cf":
            nsta = max(1, int(v["sta_len"] * self.fs))
            nlta = max(2, int(v["lta_len"] * self.fs))
            csl_full = sta_lta_centred_multichan(self.Dfull, nsta, nlta)
            self.csl_raw = csl_full[self.n_lead:, :]                     # V3
            level = "track"
        if level == "track":
            m = self._mask_vec()                                          # V3
            if self.mask_on and m is not None and m.any():                # V3
                self.csl = self.csl_raw.copy()                            # V3
                self.csl[:, m] = 0.0                                      # V3
            else:                                                         # V3
                self.csl = self.csl_raw                                   # V3
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

    # -- keys: `m` toggles the mask, everything else is v2 ------------------------------
    def _on_key(self, ev):
        from dasapp.picker import _FULLWIDTH
        k = _FULLWIDTH.get(ev.key, ev.key)
        if k == "m":
            if self.mask_full is None:
                print("[tune-v3] no channel mask loaded (run_v3.py scan first)")
                return
            self.mask_on = not self.mask_on
            self._recompute("track")
        else:
            super()._on_key(ev)

    # V3-copy of Tuner._draw — identical, plus mask stripes / status / v3 help rows.
    def _draw(self):
        ax = self.ax
        ax.clear()
        nt = self.csl.shape[0]
        t1 = nt / self.fs
        dist = (self.d1 + np.arange(self.nch) * self.dx) / 1000.0

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

        # V3: masked channels — red stripes over the full width when the mask is ON,
        # thin edge ticks when OFF (so you still see WHERE the mask would act).
        m = self._mask_vec()
        if m is not None and m.any():
            edges = np.where(np.diff(m.astype(int)) != 0)[0]
            bnds = np.concatenate([[0], edges + 1, [len(m)]])
            for i in range(len(bnds) - 1):
                lo, hi = bnds[i], bnds[i + 1] - 1
                if m[lo]:
                    if self.mask_on:
                        ax.axhspan(dist[lo], dist[hi], color="crimson",
                                   alpha=0.10, zorder=2, linewidth=0)
                    else:
                        ax.axhspan(dist[lo], dist[hi], xmax=0.008, color="crimson",
                                   alpha=0.8, zorder=2, linewidth=0)

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
        n_masked = 0 if m is None else int(m.sum())                      # V3
        mask_note = (f"mask {'ON' if self.mask_on else 'OFF'} ({n_masked} ch, key m)"
                     if n_masked else "no mask")                          # V3
        self.status = (f"tracked {len(self.wavefronts)} wavefronts · "
                       f"{npass} pass min_channels(≥{int(mc)}) · "
                       f"empty {getattr(self, 'empty_chan', 0)}/{self.nch} ch · "
                       f"{mask_note}")                                    # V3
        ax.set_title(
            f"TUNE·v3  ·  event {self.idx+1}/{len(self.cands)}  "
            f"{ev['start_time'].strftime('%Y-%m-%d %H:%M:%S')}  ·  bg: {nm}   ·   press  h  for keys\n"
            f"{self.status}   (lime = passes min_channels · grey = tracked but below it)",
            fontsize=11)
        if self.help_visible[0]:
            draw_shortcut_overlay(ax, _V3_HELP_ROWS)                      # V3
        self.fig.canvas.draw_idle()


def selftest_v3(cfg, cands, out_png, start_idx=0):
    """Headless: build the v3 tuner, toggle the mask, verify tracking reacts."""
    import os
    t = TunerV3(cfg, cands, start_idx=start_idx, interactive=False)
    n_on = len(t.wavefronts)
    class _Ev:                                        # minimal key-event stand-in
        key = "m"
    t._on_key(_Ev())
    n_off = len(t.wavefronts)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    t.fig.savefig(out_png, dpi=110)
    return {"png": out_png, "nch": t.nch, "mask_on_wavefronts": n_on,
            "mask_off_wavefronts": n_off, "status": t.status}
