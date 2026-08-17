"""Command-line entry point."""

from __future__ import annotations

import argparse
import json

from .pipeline import build


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="reasoned-mmps",
        description="Build auditable medicinal-chemistry moves with reasons.",
    )
    parser.add_argument(
        "command", choices=["build", "build-pilot"], help="Pipeline command to run"
    )
    args = parser.parse_args()
    if args.command in {"build", "build-pilot"}:
        print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
