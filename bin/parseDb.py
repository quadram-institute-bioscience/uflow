#!/usr/bin/env python3
"""
Check the database file and make a symbolic link adding the appropriate suffix
"""

import os, sys
import argparse

def dbType(file):
    """
    Check database type and return the infered extension
    """
    # Check if file is gzipped via magic number
    with open(file, 'rb') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':
        return '.fa.gz'
    else:
        return '.RData'
if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Check the database file and make a symbolic link adding the appropriate suffix")
    args.add_argument("-i", "--input", help="Input database file", required=True)
    args.add_argument("-m", "--move", help="Move the file, otherwise symlink", action="store_true")
    args = args.parse_args()

    # Output basename
    base = os.path.basename(args.input)
    ext = dbType(args.input)
    print(f"{args.input} -> {base + ext}", file=sys.stderr)
    if args.move:
        os.rename(args.input, base + ext)
    else:
        os.symlink(args.input, base + ext)