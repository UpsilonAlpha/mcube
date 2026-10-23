#!/usr/bin/env python3
"""Membrane health and partner-pose report for a membrane MD trajectory (any system).

Thin wrapper: the analysis lives in mcube.membrane_analysis (also run by `mcube postprocess` for a membrane
config). `python membrane_report.py --help` lists the options.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mcube.membrane_analysis import main  # noqa: E402

if __name__ == "__main__":
    main()
