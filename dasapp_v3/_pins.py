"""
Hash guard for the v2 files that v3 either subclasses or carries pinned copies of.

v3 duplicates a handful of v2 code blocks verbatim (documented in each module's
docstring) because v2 has no injection hook there and must not be edited. If a v2
file changes in the future, the pinned copies could silently diverge — this module
makes that loud instead: `verify()` warns (or raises, in selftest) whenever a pinned
file's SHA256 no longer matches the value recorded at v3-creation time.
"""

import os
import hashlib

_APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# recorded against v2.0.0 (git 3f1241d)
PINS = {
    "dasapp/detect.py":
        "dd523653765b590ccae4aed5908a7a12d31b32c849f2bcdb8e9e8e04b6561bc8",
    "dasapp/picker.py":
        "c24f4496446119bc8f3658985f345cc3b7d34404fc333683b2d79b82502a2a24",
    "dasapp/tune.py":
        "3ba3d7bc74793d3fee9c9ea4ac5c92a82c31ab96f930e4dcf8b936fba47f1448",
    "dasapp/config.py":
        "6631f17b0b2a64528edb84b93569b8656f8f17864f5202d9177df425057f6729",
    "vendor/functions/readers.py":
        "a194498b0eae1e8e26f713b28b6ee354c56abe9a0883e39be59d6ab1574c953d",
    "vendor/functions/stalta.py":
        "97ea3b235e2bb1162a096d674d1e33e11612daa8c648ce5ddf2fcc3868a62a7b",
}


def verify(strict=False, quiet=True):
    """Check every pinned v2 file. Returns the list of mismatched paths.

    strict=True raises RuntimeError on any mismatch (used by selftest); otherwise a
    loud warning is printed and v3 keeps running (the mismatch may be harmless, but
    the pinned copies in dasapp_v3/ must be reconciled by hand).
    """
    bad = []
    for rel, want in PINS.items():
        p = os.path.join(_APP_ROOT, rel)
        try:
            with open(p, "rb") as f:
                got = hashlib.sha256(f.read()).hexdigest()
        except OSError:
            got = "<missing>"
        if got != want:
            bad.append(rel)
    if bad:
        msg = ("[v3 pin check] v2 files changed since the v3 pinned copies were made: "
               + ", ".join(bad)
               + " — reconcile dasapp_v3/ (see V3_CHANGES.md) before trusting results.")
        if strict:
            raise RuntimeError(msg)
        print(msg)
    elif not quiet:
        print("[v3 pin check] OK — all pinned v2 files unchanged.")
    return bad
