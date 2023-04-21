#!/usr/bin/env python3
"""
Process a set of tables with two columns: IDs and values and 
produce a joined table with ID and a column for each sample
"""
import os, sys
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process a set of tables with two columns: IDs and values and produce a joined table with ID and a column for each sample"
    )
    parser.add_argument("TABS", help="Input table files", nargs="+")
    parser.add_argument("-o", "--output", help="Output table file", required=True)
    parser.add_argument("-s", "--sample-from-header", help="Sample names from table headers", action="store_true")
    parser.add_argument("--sort", help="Sort sample names", action="store_true")
    parser.add_argument("--otuid", help="OTU ID column name", default="#OTUID")
    parser.add_argument(
        "-v", "--verbose", help="Verbose output", action="store_true"
    )
    args = parser.parse_args()

    if args.verbose:
        print(args)
    
    # Read input files
    otus = []
    samples = []
    tables = {}
    for singleTabFile in args.TABS:
        table = {}
        
        samplename = os.path.basename(singleTabFile).split(".")[0]
        
        if args.verbose:
            print(f"Reading {singleTabFile}")
        with open(singleTabFile, "r") as f:
            c=0
            for line in f:
                c+=1
                if c==1:
                    if args.sample_from_header:
                        samplename = line.strip().split("\t")[1]
                    continue
                line = line.strip()
                if line:
                    id, val = line.split()
                    if not id in otus:
                        otus.append(id)
                    table[id] = val
        tables[samplename] = table 
        samples.append(samplename)

    if args.sort:
        samples.sort()
    # Join tables
    with open(args.output, "w") as f:
        # Header
        f.write(f"{args.otuid}\t")
        for samplename in samples:
            sep = "\t" if samplename != samples[-1] else ""
            f.write(f"{samplename}{sep}")
        f.write("\n")

        print(f"Last sample: {samples[-1]}", file=sys.stderr)
        for id in otus:
            f.write(f"{id}\t")
            for samplename in samples:
                sep = "\t"
                if samplename == samples[-1]:
                    sep = "\n"
                if id in tables[samplename]:
                    # Check if id is the last one

                    
                    f.write(f"{tables[samplename][id]}{sep}")
                else:
                    f.write(f"0{sep}")
                
            #f.write("\n")
