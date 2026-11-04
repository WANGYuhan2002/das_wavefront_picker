#!/bin/bash
# Double-click this file in Finder to launch das_wavefront_picker.
# It opens a Terminal, runs a bundled-data demo (detect + interactive picker).
# For your own data, edit config.json and run:  python3 run.py both
cd "$(dirname "$0")"
echo "=== das_wavefront_picker ==="
echo "Launching the bundled demo (detect on the sample, then the interactive picker)…"
echo "For your own data: edit config.json, then run  python3 run.py both"
echo
exec /usr/bin/env python3 run.py demo "$@"
