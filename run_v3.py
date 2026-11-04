#!/usr/bin/env python3
"""
das_wavefront_picker v3 — the v2 pipeline + a channel-quality (bad-channel) workflow.
v2 (`run.py`) is completely untouched and keeps working; v3 is additive. See V3_CHANGES.md.

    python run_v3.py scan              # channel-quality scan -> mask + report figure
    python run_v3.py mask              # re-derive the mask from a scan (new thresholds, no re-read)
    python run_v3.py detect            # masked STA/LTA detection -> Output/v3_<exp>_<mask>/
    python run_v3.py pick              # v3 picker: mask stripes + event-local SNR conf gate
    python run_v3.py tune              # v3 tuner: `m` toggles the mask live
    python run_v3.py both              # detect, then pick
    python run_v3.py demo              # bundled sample: scan + detect + pick  (--tune / --quick)
    python run_v3.py info              # resolved config + mask + catalog stats
    python run_v3.py selftest          # headless end-to-end check (scan/detect/pick/tune)

Typical first session on your own data:
    python run_v3.py scan  --data /archive/DAS/ --start ... --end ...
      -> open Output/ChannelQuality/scan_*/quality_report.png, eyeball the mask,
         hand-edit channel_mask.csv if needed (or set quality.manual_bad/manual_good)
    python run_v3.py detect            # uses the newest compatible mask automatically
    python run_v3.py pick
"""

import os
import sys
import argparse

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from dasapp.config import add_common_cli, DEMO_DATA_DIR              # v2, untouched
from dasapp_v3.config_v3 import add_v3_cli, config_v3_from_args      # noqa: E402
from dasapp_v3 import _pins                                          # noqa: E402


def _cfg(args, need_data=True):
    cfg = config_v3_from_args(args)
    if need_data and not cfg.das_data_dir:
        sys.exit("[error] no DAS data directory set.\n"
                 "        Set \"das_data_dir\" in config.json (or pass --data <dir>).")
    return cfg


def _resolve_mask(cfg, args, require=False):
    from dasapp_v3.quality import resolve_mask
    resolve_mask(cfg, no_mask=getattr(args, "no_mask", False), require=require)
    if cfg.chan_mask is not None:
        print(f">>> channel mask: {int(cfg.chan_mask.sum())} ch masked "
              f"({cfg.chan_mask_path})")
    else:
        print(">>> channel mask: none (v2-equivalent detection)")
    return cfg


# ------------------------------------------------------------------ subcommands -------
def _cmd_scan(args):
    from dasapp_v3.quality import run_scan
    cfg = _cfg(args)
    print(cfg.summary(), "\n")
    scan_dir, stats = run_scan(cfg)
    print("\nnext: inspect quality_report.png / channel_quality.csv, hand-edit "
          "channel_mask.csv if needed, then `python run_v3.py detect`.")
    return scan_dir


def _cmd_mask(args):
    """Re-derive flags+mask from an existing scan with (possibly new) thresholds."""
    import glob
    from dasapp_v3.quality import rethreshold, rewrite_artifacts
    cfg = _cfg(args, need_data=False)
    scan_dir = args.scan_dir
    if not scan_dir:
        cands = sorted(glob.glob(os.path.join(cfg.quality_root, "scan_*")),
                       key=os.path.getmtime, reverse=True)
        if not cands:
            sys.exit(f"[error] no scans under {cfg.quality_root} — run `run_v3.py scan`.")
        scan_dir = cands[0]
    print(f">>> re-thresholding {scan_dir}")
    stats, dist_km, E = rethreshold(scan_dir, cfg.quality)
    rewrite_artifacts(scan_dir, cfg.quality, stats, dist_km, E)
    from dasapp_v3.quality import print_scan_summary
    print_scan_summary(cfg, stats, len(dist_km))


def _cmd_detect(args):
    from dasapp_v3.detect_v3 import run_detection_v3
    from dasapp_v3.quality import mask_guardrails
    cfg = _cfg(args)
    _resolve_mask(cfg, args, require=not getattr(args, "no_mask", False))
    if cfg.chan_mask is not None:
        for w in mask_guardrails(cfg.chan_mask, cfg):
            print(f">>> WARNING: {w}")
    print(cfg.summary(), "\n")
    path = run_detection_v3(cfg, resume=not args.no_resume)
    return cfg, path


def _cmd_compare(args):
    from dasapp_v3.compare_ab import run_compare
    cfg = _cfg(args, need_data=False)
    run_compare(args.a, args.b, tol_s=float(args.tol), fs=cfg.fs)


def _cmd_pick(args, cfg=None, catalog=None):
    from dasapp import catalog as cat                                 # v2, untouched
    from dasapp_v3.picker_v3 import PickerV3
    if cfg is None:
        cfg = _cfg(args)
        _resolve_mask(cfg, args, require=False)
    catalog = catalog or getattr(args, "catalog", None) or cfg.catalog_path
    if not os.path.exists(catalog):
        sys.exit(f"[error] catalog not found: {catalog}\n"
                 f"        run `python run_v3.py detect` first, or pass --catalog <path>.")
    exp_dir = os.path.dirname(os.path.abspath(catalog))
    events, schema = cat.load_catalog(catalog)
    print(f">>> {len(events)} candidate wavefronts ({schema}) from {catalog}")
    if not events:
        sys.exit("[error] catalog is empty — nothing to pick.")
    start = int(getattr(args, "event", 0) or 0)
    PickerV3(cfg, events, exp_dir, start_idx=start, interactive=True).run()


def _cmd_both(args):
    cfg, path = _cmd_detect(args)
    _cmd_pick(args, cfg=cfg, catalog=path)


def _cmd_tune(args):
    from dasapp import catalog as cat
    from dasapp_v3.tune_v3 import TunerV3
    cfg = _cfg(args)
    _resolve_mask(cfg, args, require=False)
    catalog = getattr(args, "catalog", None) or cfg.catalog_path
    if not os.path.exists(catalog):
        sys.exit(f"[error] need a catalog to pick preview windows from: {catalog}\n"
                 f"        run `python run_v3.py detect` first, or pass --catalog <path>.")
    cands, schema = cat.load_catalog(catalog)
    if not cands:
        sys.exit("[error] catalog is empty.")
    print(f">>> {len(cands)} candidate windows ({schema}) from {catalog}")
    TunerV3(cfg, cands, start_idx=int(getattr(args, "event", 0) or 0),
            interactive=True).run()


def _cmd_demo(args):
    if not DEMO_DATA_DIR:
        sys.exit("[error] demo data not found (see README §1).")
    if not getattr(args, "das_data_dir", None):
        args.das_data_dir = DEMO_DATA_DIR
    args.method = args.method or "B"
    args.t_start_str = args.t_start_str or "2024-04-25T15:00:00"
    args.t_end_str = args.t_end_str or ("2024-04-25T15:05:00" if args.quick
                                        else "2024-04-25T15:30:00")
    # scan first (unless a mask already exists / --no-mask)
    from dasapp_v3.quality import find_latest_mask, run_scan
    cfg = _cfg(args)
    if not getattr(args, "no_mask", False) and find_latest_mask(cfg) is None:
        print(">>> demo: no channel mask yet — running the quality scan first\n")
        run_scan(cfg)
    cfg, path = _cmd_detect(args)
    if args.no_pick:
        print(f"\n[demo] detection done. Pick: python run_v3.py pick --catalog {path}")
        return
    if getattr(args, "tune", False):
        from dasapp import catalog as cat
        from dasapp_v3.tune_v3 import TunerV3
        cands, _schema = cat.load_catalog(path)
        TunerV3(cfg, cands, start_idx=int(getattr(args, "event", 0) or 0),
                interactive=True).run()
    else:
        _cmd_pick(args, cfg=cfg, catalog=path)


def _cmd_info(args):
    from dasapp import catalog as cat
    cfg = _cfg(args, need_data=False)
    _resolve_mask(cfg, args, require=False)
    print(cfg.summary())
    if cfg.chan_mask is not None:
        ch = list(map(int, cfg.chan_mask.nonzero()[0]))
        print(f"mask: {len(ch)} channels -> {ch}")
    cpath = getattr(args, "catalog", None) or cfg.catalog_path
    print(f"\ncatalog: {cpath}")
    if os.path.exists(cpath):
        events, schema = cat.load_catalog(cpath)
        src = {}
        for e in events:
            src[e["source"]] = src.get(e["source"], 0) + 1
        print(f"  {len(events)} events (schema {schema}); by method: {src}")
    else:
        print("  (not created yet — run detect)")


def _cmd_selftest(args):
    """Headless end-to-end: scan -> masked detect -> v3 pick export -> v3 tune toggle."""
    from dasapp import catalog as cat
    from dasapp_v3.quality import run_scan, resolve_mask
    from dasapp_v3.detect_v3 import run_detection_v3
    from dasapp_v3.picker_v3 import selftest_v3 as pick_selftest
    from dasapp_v3.tune_v3 import selftest_v3 as tune_selftest
    if DEMO_DATA_DIR and not args.das_data_dir:
        args.das_data_dir = DEMO_DATA_DIR
        args.t_start_str = args.t_start_str or "2024-04-25T15:00:00"
        args.t_end_str = args.t_end_str or "2024-04-25T15:10:00"
        args.method = args.method or "B"
    cfg = _cfg(args)
    _pins.verify(strict=True)
    print("[selftest] pin check OK (v2 files unchanged)")
    scan_dir, stats = run_scan(cfg, progress=False)
    print(f"[selftest] scan OK: {int(stats['mask'].sum())} channels masked -> {scan_dir}")
    resolve_mask(cfg, no_mask=False, require=True)
    if not os.path.exists(cfg.catalog_path):
        run_detection_v3(cfg, progress=False, resume=False)
    events, schema = cat.load_catalog(cfg.catalog_path)
    if not events:
        sys.exit("[selftest] no events detected in the test window.")
    out_png = os.path.join(cfg.exp_dir, "selftest_v3_wavefront.png")
    info = pick_selftest(cfg, events, cfg.exp_dir, out_png, idx=0)
    print(f"[selftest] pick OK  candidates={len(events)} schema={schema}")
    for k, v in info.items():
        print(f"   {k}: {v}")
    tinfo = tune_selftest(cfg, events, os.path.join(cfg.exp_dir, "selftest_v3_tune.png"))
    print(f"[selftest] tune OK  mask_on={tinfo['mask_on_wavefronts']} "
          f"mask_off={tinfo['mask_off_wavefronts']} wavefronts")
    print("[selftest] OK")


def main():
    ap = argparse.ArgumentParser(prog="run_v3.py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("scan", help="channel-quality scan -> mask + report")
    add_common_cli(p); add_v3_cli(p)
    p.set_defaults(func=_cmd_scan)

    p = sub.add_parser("mask", help="re-derive the mask from an existing scan")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--scan-dir", help="scan folder (default: newest under ChannelQuality/)")
    p.set_defaults(func=_cmd_mask)

    p = sub.add_parser("detect", help="masked STA/LTA detection")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--no-resume", action="store_true")
    p.set_defaults(func=_cmd_detect)

    p = sub.add_parser("pick", help="v3 interactive picker")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--catalog"); p.add_argument("--event", type=int, default=0)
    p.set_defaults(func=_cmd_pick)

    p = sub.add_parser("both", help="detect then pick")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--event", type=int, default=0); p.add_argument("--catalog")
    p.set_defaults(func=_cmd_both)

    p = sub.add_parser("tune", help="v3 parameter tuner (mask toggle: m)")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--catalog"); p.add_argument("--event", type=int, default=0)
    p.set_defaults(func=_cmd_tune)

    p = sub.add_parser("demo", help="bundled sample: scan + detect + pick")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--quick", action="store_true", help="5-min window instead of 30")
    p.add_argument("--tune", action="store_true", help="open the tuner instead of the picker")
    p.add_argument("--no-pick", action="store_true", help="detect only")
    p.add_argument("--event", type=int, default=0); p.add_argument("--catalog")
    p.set_defaults(func=_cmd_demo)

    p = sub.add_parser("compare", help="A/B-compare two catalogs (e.g. v2 vs v3)")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--a", required=True, help="baseline catalog.csv (e.g. v2)")
    p.add_argument("--b", required=True, help="comparison catalog.csv (e.g. v3)")
    p.add_argument("--tol", default=1.0, help="event match tolerance (s)")
    p.set_defaults(func=_cmd_compare)

    p = sub.add_parser("info", help="config + mask + catalog stats")
    add_common_cli(p); add_v3_cli(p)
    p.add_argument("--catalog")
    p.set_defaults(func=_cmd_info)

    p = sub.add_parser("selftest", help="headless end-to-end check")
    add_common_cli(p); add_v3_cli(p)
    p.set_defaults(func=_cmd_selftest)

    args = ap.parse_args()
    if args.cmd != "selftest":
        _pins.verify(strict=False)          # warn loudly if v2 files drifted
    args.func(args)


if __name__ == "__main__":
    main()
