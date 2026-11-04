#!/usr/bin/env python3
"""
das_wavefront_picker — one command for STA/LTA detection + interactive wavefront picking.

    python run.py demo                 # bundled-sample demo of the PICKER  (detect + pick)
    python run.py demo --tune          # bundled-sample demo of the parameter TUNER (detect + tune)
    python run.py detect               # STA/LTA detection -> Output/<exp>/catalog.csv
    python run.py pick                 # open the interactive per-event picker on that catalog
    python run.py tune                 # open the detection-parameter tuner (sliders) on that catalog
    python run.py both                 # detect, then immediately pick
    python run.py info                 # show the resolved config + catalog stats
    python run.py selftest             # headless end-to-end check (no window; renders a PNG)

Point it at your own data by editing config.json (das_data_dir, fs, dx, band, range, time
range, method, thresholds) or with flags, e.g.:

    python run.py detect --data /archive/DATA_DEC_80Hz_16m_SORTED/ --fs 80 --dx 16.335 \
                         --band 5 20 --range 8000 14000 --method B \
                         --start 2024-04-25T15:00:00 --end 2024-04-25T15:30:00

Open an EXISTING matched_filter_release catalog directly (method A v1 or method B v2):

    python run.py pick --catalog /path/to/sta_lta_catalog_v2.csv
"""

import os
import sys
import argparse

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from dasapp.config import Config, add_common_cli, config_from_args, DEMO_DATA_DIR  # noqa: E402


def _cmd_detect(args):
    from dasapp.detect import run_detection
    cfg = config_from_args(args)
    if not cfg.das_data_dir:
        sys.exit("[error] no DAS data directory set.\n"
                 "        Set \"das_data_dir\" in config.json (or pass --data <dir>) to a folder\n"
                 "        of DAS .h5 files. To run the built-in demo, drop the sample h5 files\n"
                 "        into das_wavefront_picker/demo_data/ (see README section 1).")
    print(cfg.summary(), "\n")
    path = run_detection(cfg, resume=not args.no_resume)
    return cfg, path


def _cmd_pick(args, cfg=None, catalog=None):
    from dasapp import catalog as cat
    from dasapp.picker import Picker
    if cfg is None:
        cfg = config_from_args(args)
    catalog = catalog or getattr(args, "catalog", None) or cfg.catalog_path
    if not os.path.exists(catalog):
        sys.exit(f"[error] catalog not found: {catalog}\n"
                 f"        run `python run.py detect` first, or pass --catalog <path>.")
    exp_dir = os.path.dirname(os.path.abspath(catalog))
    events, schema = cat.load_catalog(catalog)
    print(f">>> {len(events)} candidate wavefronts ({schema}) from {catalog}")
    if not events:
        sys.exit("[error] catalog is empty — nothing to pick.")
    start = int(getattr(args, "event", 0) or 0)
    Picker(cfg, events, exp_dir, start_idx=start, interactive=True).run()


def _cmd_both(args):
    cfg, path = _cmd_detect(args)
    _cmd_pick(args, cfg=cfg, catalog=path)


def _cmd_demo(args):
    if not DEMO_DATA_DIR:
        sys.exit("[error] demo data not found.\n"
                 "        Put the sample DAS .h5 files into das_wavefront_picker/demo_data/\n"
                 "        (see README section 1), or set das_data_dir in config.json and run\n"
                 "        `python run.py detect` on your own data.")
    # demo defaults — but honour an explicit --data / --method / window if the user gave one
    if not getattr(args, "das_data_dir", None):
        args.das_data_dir = DEMO_DATA_DIR
    args.method = args.method or "B"
    args.t_start_str = args.t_start_str or "2024-04-25T15:00:00"
    args.t_end_str = args.t_end_str or ("2024-04-25T15:05:00" if args.quick
                                        else "2024-04-25T15:30:00")
    cfg, path = _cmd_detect(args)
    if args.no_pick:
        print(f"\n[demo] detection done. Pick: python run.py pick --catalog {path}\n"
              f"                       Tune: python run.py tune --catalog {path}")
        return
    if getattr(args, "tune", False):                    # demo the parameter tuner instead of the picker
        from dasapp import catalog as cat
        from dasapp.tune import Tuner
        cands, _schema = cat.load_catalog(path)
        print(">>> demo: opening the detection-parameter TUNER (sliders). "
              "Save to config, then run.py detect.")
        Tuner(cfg, cands, start_idx=int(getattr(args, "event", 0) or 0), interactive=True).run()
    else:
        _cmd_pick(args, cfg=cfg, catalog=path)


def _cmd_tune(args):
    """Interactive detection-parameter tuner: sliders re-run STA/LTA + tracking live on a
    catalog window; Save writes params back to config.json for the full `detect` batch."""
    from dasapp import catalog as cat
    from dasapp.tune import Tuner
    cfg = config_from_args(args)
    if not cfg.das_data_dir:
        sys.exit("[error] set das_data_dir (config.json) or pass --data — tuning needs the raw DAS.")
    catalog = getattr(args, "catalog", None) or cfg.catalog_path
    if not os.path.exists(catalog):
        sys.exit(f"[error] need a catalog to pick preview windows from: {catalog}\n"
                 f"        run `python run.py detect` first, or pass --catalog <path>.")
    cands, schema = cat.load_catalog(catalog)
    if not cands:
        sys.exit("[error] catalog is empty.")
    print(f">>> {len(cands)} candidate windows ({schema}) from {catalog}")
    Tuner(cfg, cands, start_idx=int(getattr(args, "event", 0) or 0), interactive=True).run()


def _cmd_info(args):
    from dasapp import catalog as cat
    cfg = config_from_args(args)
    print(cfg.summary())
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
    """Headless: ensure a small catalog exists (detect a short window if needed),
    then run the picker's non-interactive selftest on event 0."""
    from dasapp import catalog as cat
    from dasapp.detect import run_detection
    from dasapp.picker import selftest
    if DEMO_DATA_DIR and not args.das_data_dir:
        args.das_data_dir = DEMO_DATA_DIR
        args.t_start_str = args.t_start_str or "2024-04-25T15:00:00"
        args.t_end_str = args.t_end_str or "2024-04-25T15:05:00"
        args.method = args.method or "B"
    cfg = config_from_args(args)
    if not os.path.exists(cfg.catalog_path):
        run_detection(cfg, resume=False)
    events, schema = cat.load_catalog(cfg.catalog_path)
    if not events:
        sys.exit("[selftest] no events detected in the test window.")
    out_png = os.path.join(cfg.exp_dir, "selftest_wavefront.png")
    info = selftest(cfg, events, cfg.exp_dir, out_png, idx=0)
    print(f"[selftest] OK  candidates={len(events)} schema={schema}")
    for k, v in info.items():
        print(f"   {k}: {v}")


def main():
    ap = argparse.ArgumentParser(prog="run.py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("detect", help="run STA/LTA detection")
    add_common_cli(p); p.add_argument("--no-resume", action="store_true")
    p.set_defaults(func=_cmd_detect)

    p = sub.add_parser("pick", help="interactive per-event picker")
    add_common_cli(p); p.add_argument("--catalog"); p.add_argument("--event", type=int, default=0)
    p.set_defaults(func=_cmd_pick)

    p = sub.add_parser("both", help="detect then pick")
    add_common_cli(p); p.add_argument("--no-resume", action="store_true")
    p.add_argument("--event", type=int, default=0); p.add_argument("--catalog")
    p.set_defaults(func=_cmd_both)

    p = sub.add_parser("demo", help="run on the bundled sample (detect + open the PICKER)")
    add_common_cli(p); p.add_argument("--no-resume", action="store_true")
    p.add_argument("--quick", action="store_true", help="5-min window instead of 30")
    p.add_argument("--tune", action="store_true", help="open the parameter TUNER instead of the picker")
    p.add_argument("--no-pick", action="store_true", help="detect only, don't open any window")
    p.add_argument("--event", type=int, default=0); p.add_argument("--catalog")
    p.set_defaults(func=_cmd_demo)

    p = sub.add_parser("tune", help="interactive detection-parameter tuner (sliders)")
    add_common_cli(p); p.add_argument("--catalog"); p.add_argument("--event", type=int, default=0)
    p.set_defaults(func=_cmd_tune)

    p = sub.add_parser("info", help="show config + catalog stats")
    add_common_cli(p); p.add_argument("--catalog")
    p.set_defaults(func=_cmd_info)

    p = sub.add_parser("selftest", help="headless end-to-end check")
    add_common_cli(p); p.set_defaults(func=_cmd_selftest)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
