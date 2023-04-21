#!/usr/bin/env python
"""
Prepare the FASTQ reads for USEARCH
August 2022, Andrea Telatin
"""

import os, sys
import argparse

def main():
    args = argparse.ArgumentParser("Prepare the FASTQ reads for USEARCH")
    args.add_argument("FASTQ", help="FASTQ files", nargs="+")
    args.add_argument("-o", "--output", help="Output directory [default: %(default)s]", default="usearch-reads")
if __name__=="__main__":
    main()