#!/usr/bin/env python3
"""Run the frozen reasoned-MMP pilot."""

from __future__ import annotations

import json

from reasoned_mmp.pipeline import build


if __name__ == "__main__":
    print(json.dumps(build(), indent=2, sort_keys=True))
