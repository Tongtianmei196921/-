#!/usr/bin/env python3
"""Compatibility CLI that forwards to the packaged implementation."""

from .drugreflector.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
