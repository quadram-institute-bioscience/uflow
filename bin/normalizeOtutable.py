#!/usr/bin/env python3
import os, sys
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Normalize OTU table')
    parser.add_argument('-i', '--input', help='Input OTU table', required=True)
    parser.add_argument('-o', '--output', help='Output OTU table', required=True)
    parser.add_argument('-s', '--separator', help='Separator', default='\t')
     
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print('Error: OTU table file not found: ' + args.input)
        sys.exit(1)
    
    output_file = args.output
    
    df = pd.read_csv(args.input, sep=args.separator, index_col=0)
    
    # Normalize by total of each column
    df_norm = df.div(df.sum(axis=0), axis=1)

    
    df_norm.to_csv(output_file, sep='\t')