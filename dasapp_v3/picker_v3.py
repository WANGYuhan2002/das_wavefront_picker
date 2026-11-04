"""
v3 interactive picker = the v2 Picker with channel quality wired into per-channel
confidence, without modifying dasapp/picker.py.

Per-channel confidence rule (the heart of it):

    conf = static_ok  AND  edge_ok  AND  ( (linked AND consistent)  OR  cf_gate )

  static_ok : not masked by the quality scan (masked channels tint red in the window,
              distinct from the ordinary grey low-confidence bands);
  edge_ok   : honours cfg.edge_margin_km (declared in v2 config but never consumed);
  linked    : the channel is an actual member of the method-B tracked chain — a linked
              channel passes UNCONDITIONALLY, because the tracker's moveout-consistency
              is itself evidence (the low-SNR-but-consistent escape hatch);
  consistent: linked members must not be gross moveout outliers (mislink prune):
              |arrival − local median| <= max(prune_k·1.4826·MAD, 3 samples) over a
              prune_win_ch moving window;
  cf_gate   : unlinked in-span channels are re-admitted if the local CF peak within
              ±event_snr_win_s of the seed arrival exceeds event_snr_gate·threshold_level
              — v2's own method-A convention (picker._seed_auto), applied to method B.

Why: v2 marks EVERY channel inside a tracked span conf=True, but a median 29% of
in-span channels were never linked (measured on the demo baseline), and on weak events
40-60% of those fail any reasonable SNR test — they simply don't see the event (the
user's second notion of channel SNR). Validated: 100% of linked channels also pass the
cf_gate on the demo (the gate never fights the tracker), while on weak events it
correctly rejects the blind in-span channels.

The operator keeps the last word: `c` still toggles confidence bands, and every manual
change is recorded as qc_code=3 in the exported h5.

Export additions (after `super().save()`, appended via h5py — the v2/release schema is
untouched, `conf` keeps its boolean contract): per-channel `qc_code` (0 good ·
1 static-masked · 2 event-local reject · 3 manually toggled) and `snr_db` (post/pre
energy ratio around the final pick; diagnostic/weighting metadata, never a gate).

Only two v2 methods are copied (marked `# V3-copy`): Picker.load (to construct
WaveWindowV3 — v2 has no factory hook) and the wiggle drawer (to tint masked traces).
"""

import os
import copy
import numpy as np
from scipy.ndimage import median_filter

from dasapp.picker import Picker, WaveWindow            # v2, untouched

try:
    import h5py
except Exception:
    h5py = None


# ============================================================= window =================
class WaveWindowV3(WaveWindow):
    """WaveWindow + the v3 confidence rule applied to the seed."""

    def __init__(self, cand, cfg, exp_dir):
        self._linked = None                     # bool[nch] chain membership (method B)
        self.snr_ev = None                      # per-channel local CF peak
        super().__init__(cand, cfg, exp_dir)    # builds D/CF_raw and calls _seed_*
        self._apply_quality_gates()

    # -- seed hooks: record chain membership, then defer to v2 --------------------------
    def _seed_from_tracked(self, rel):
        seed, conf = super()._seed_from_tracked(rel)
        ch = np.clip(rel[0].astype(int), 0, self.nch - 1)
        self._linked = np.zeros(self.nch, bool)
        self._linked[ch] = True
        return seed, conf

    def _seed_auto(self):
        self._linked = None                     # method A: no chain, keep v2's val gate
        return super()._seed_auto()

    # -- the v3 confidence rule ---------------------------------------------------------
    def _apply_quality_gates(self):
        cfg = self.cfg
        q = getattr(cfg, "quality", None) or {}
        gate = float(q.get("event_snr_gate", 0.5)) * float(cfg.threshold_level)
        win = float(q.get("event_snr_win_s", 0.5))

        # event-local CF peak within ±win of the (dense) seed arrival, per channel
        w = max(1, int(win * self.fs))
        idx = np.clip((self.seed * self.fs).astype(int), 0, self.nt - 1)
        snr = np.empty(self.nch)
        for k in range(self.nch):
            i0, i1 = max(0, idx[k] - w), min(self.nt, idx[k] + w + 1)
            col = self.CF_raw[i0:i1, k]
            snr[k] = float(np.nanmax(col)) if col.size else 0.0
        self.snr_ev = snr
        cf_gate = snr > gate

        if self._linked is not None:            # method B: (linked & consistent) | cf_gate
            consistent = self._prune_consistent(q)
            self.seed_conf = self.seed_conf & ((self._linked & consistent) | cf_gate)
        else:                                    # method A: v2's own val gate stands
            self.seed_conf = self.seed_conf & cf_gate

        # static mask: quality-scan channels can never be confident
        mask = getattr(cfg, "chan_mask", None)
        if mask is not None and len(mask) >= self.nch:
            self.mask = np.asarray(mask[: self.nch], bool)
            self.seed_conf = self.seed_conf & ~self.mask
        else:
            self.mask = np.zeros(self.nch, bool)

        # edge margin (v2 config key, consumed here)
        e = float(getattr(cfg, "edge_margin_km", 0.0) or 0.0)
        if e > 0:
            d = self.dist_km
            self.seed_conf = self.seed_conf & (d >= d.min() + e) & (d <= d.max() - e)

    def _prune_consistent(self, q):
        """Mislink prune on linked members: True where the arrival is locally consistent."""
        consistent = np.ones(self.nch, bool)
        k = int(q.get("prune_win_ch", 15))
        if k <= 0 or self._linked.sum() < 5:
            return consistent
        members = np.where(self._linked)[0]
        arr = self.seed[members]                 # display seconds at the linked channels
        size = min(k, len(arr))
        if size % 2 == 0:
            size -= 1
        if size < 3:
            return consistent
        med = median_filter(arr, size=size, mode="nearest")
        res = arr - med
        madv = np.median(np.abs(res - np.median(res)))
        tol = max(float(q.get("prune_k", 3.0)) * 1.4826 * madv, 3.0 / self.fs)
        consistent[members[np.abs(res) > tol]] = False
        return consistent


# ============================================================= picker =================
class PickerV3(Picker):

    # V3-copy of Picker.load (dasapp/picker.py) — identical except WaveWindow ->
    # WaveWindowV3 (v2 constructs the window class inline; no smaller seam exists)
    # and one added line stashing the auto conf for manual-toggle bookkeeping.
    def load(self, idx):
        idx = int(np.clip(idx, 0, len(self.cands) - 1))
        if idx not in self._cache:
            if len(self._cache) > 4:
                self._cache.pop(next(iter(self._cache)))
            self._cache[idx] = WaveWindowV3(self.cands[idx], self.cfg, self.exp_dir)  # V3
        self.idx = idx
        self.win = self._cache[idx]
        dist = self.win.dist_km

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
        self._conf_auto = conf0.copy()                    # V3: for qc_code bookkeeping
        self.dirty = False
        self.redraw()

    # -- draw: v2 redraw + quality overlays --------------------------------------------
    def redraw(self):
        super().redraw()
        w = self.win
        mask = getattr(w, "mask", None)
        if mask is None or not mask.any():
            return
        ax = self.ax
        # red stripes: statically masked channels (distinct from grey low-conf bands)
        edges = np.where(np.diff(mask.astype(int)) != 0)[0]
        bnds = np.concatenate([[0], edges + 1, [len(mask)]])
        for i in range(len(bnds) - 1):
            lo, hi = bnds[i], bnds[i + 1] - 1
            if mask[lo]:
                ax.axhspan(w.dist_km[lo], w.dist_km[hi], color="crimson",
                           alpha=0.10, zorder=1.5, linewidth=0)
        ax.text(0.008, 0.012, f"quality mask: {int(mask.sum())}/{w.nch} ch (red)",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8.5, color="crimson", alpha=0.85, zorder=6)
        self.fig.canvas.draw_idle()

    # V3-copy of Picker._draw_wiggle_bg — identical except masked traces draw in red
    # so a bad channel's waveform is visibly flagged instead of silently blending in.
    def _draw_wiggle_bg(self, ax, gmul=1.0):
        w = self.win
        mask = getattr(w, "mask", np.zeros(w.nch, bool))
        k0, k1 = self._wiggle_krange()
        n = min(self.n_wig, k1 - k0 + 1)
        ksel = np.unique(np.linspace(k0, k1, max(n, 1)).astype(int))
        d = w.dist_km[ksel]
        if self.focus[0] and self.support is not None:
            vlo, vhi = self.support
        else:
            vlo, vhi = w.dist_km.min(), w.dist_km.max()
        n_show = int(np.sum((d >= vlo) & (d <= vhi)))
        self._n_shown = n_show
        gain = 0.45 * (vhi - vlo) / max(n_show, 1) * gmul
        for ki, dk in zip(ksel, d):
            tr = w.D[:, ki].astype(float)
            tr = tr / (np.max(np.abs(tr)) + 1e-9) * gain
            col = "crimson" if mask[ki] else "k"                       # V3
            alpha = 0.55 if mask[ki] else 1.0                          # V3
            ax.plot(w.t, dk + tr, color=col, lw=0.4, zorder=0.5, alpha=alpha)
            ax.fill_between(w.t, dk, dk + tr, where=(tr > 0), color=col,
                            linewidth=0, zorder=0.5, alpha=alpha)

    # -- export: v2 save + appended qc datasets -----------------------------------------
    def save(self):
        fname = super().save()
        if h5py is None:
            return fname
        w = self.win
        cfg = self.cfg
        q = getattr(cfg, "quality", None) or {}
        try:
            wf = self.wavefront()                           # display seconds per channel
            n0, n1 = q.get("snr_noise_s", [-2.5, -0.5])
            s0, s1 = q.get("snr_sig_s", [0.0, 0.5])
            snr_db = np.full(w.nch, np.nan, dtype=np.float32)
            for k in range(w.nch):
                t = wf[k]
                ia = (max(0, int((t + n0) * w.fs)), max(0, int((t + n1) * w.fs)))
                ib = (max(0, int((t + s0) * w.fs)), min(w.nt, int((t + s1) * w.fs)))
                if ia[1] - ia[0] >= 4 and ib[1] - ib[0] >= 4:
                    pn = float(np.mean(w.D[ia[0]:ia[1], k] ** 2))
                    ps = float(np.mean(w.D[ib[0]:ib[1], k] ** 2))
                    if pn > 0:
                        snr_db[k] = 10.0 * np.log10(ps / pn + 1e-30)
            mask = getattr(w, "mask", np.zeros(w.nch, bool))
            qc = np.zeros(w.nch, dtype=np.uint8)            # 0 good
            qc[~self.conf & ~mask] = 2                      # event-local reject
            qc[mask] = 1                                    # static mask
            qc[self.conf != self._conf_auto] = 3            # operator override
            path = os.path.join(self.picks_dir, fname)
            with h5py.File(path, "a") as f:
                for name, arr in (("qc_code", qc), ("snr_db", snr_db)):
                    if name in f:
                        del f[name]
                    f.create_dataset(name, data=arr)
                f.attrs["n_static_masked"] = int(mask.sum())
                f.attrs["n_event_low"] = int((qc == 2).sum())
                f.attrs["qc_mask_hash"] = getattr(cfg, "mask_tag", "nomask")
                f.attrs["qc_artifact"] = getattr(cfg, "chan_mask_path", "")
                f.attrs["event_snr_gate"] = float(q.get("event_snr_gate", 0.5))
                f.attrs["event_snr_win_s"] = float(q.get("event_snr_win_s", 0.5))
        except Exception as e:                              # qc metadata is best-effort
            print(f"[pick-v3] qc metadata append failed: {e}")
        return fname


# ============================================================= headless selftest ======
def selftest_v3(cfg, cands, exp_dir, out_png, idx=0):
    """Non-interactive smoke test of the v3 picker (mirrors the v2 selftest)."""
    p = PickerV3(cfg, cands, exp_dir, start_idx=idx, interactive=False)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    p.fig.savefig(out_png, dpi=120)
    fname = p.save()
    w = p.win
    return {"png": out_png, "file": fname, "label": p.label[0],
            "nt": w.nt, "nch": w.nch,
            "n_conf": int(p.conf.sum()),
            "n_masked": int(getattr(w, "mask", np.zeros(1)).sum()),
            "n_linked": int(w._linked.sum()) if w._linked is not None else -1,
            "cf_max": float(np.nanmax(w.CF_raw))}
