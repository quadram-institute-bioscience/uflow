#!/usr/bin/env python3
"""
Read an OTU table and a taxonomy file, then add a "Taxonomy" column to the OTU table
"""
#Zotu1   Root [rootrank, 100.0%]; Fungi [kingdom, 100.0%]; Ascomycota [phylum, 100.0%]; Saccharomycetes [class, 100.0%]; Saccharomycetales [order, 100.0%]; Debaryomycetaceae [family, 100.0%]; Debaryomyces [genus, 98.8%]; unclassified_Debaryomyces [species, 98.8%]
import pandas as pd
import os, sys
import argparse

def loadTaxonomy(taxonomy, fasta):
    # Load the OTU names from the fasta file
    otu_names = []
    with open(fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                otu_names.append(line.strip().split(" ")[0][1:])
    # Load the taxonomy file
    tax = pd.read_csv(taxonomy, sep=" ", index_col=0)
    
    # Check that the index has the size of the OTU names
    if len(otu_names) != len(tax.index):
        print("Error: OTU names and taxonomy file do not have the same number of OTUs")
        sys.exit(1)
    # Rename the index to the OTU names
    tax.index = otu_names
    return tax

def condenseTaxonomy(tax):
    # Condense the taxonomy to a string like k__Kingdom; p__Phylum; c__Class; o__Order; f__Family; g__Genus; s__Species
    tax = tax.copy()
    tax["Taxonomy"] = tax.apply(lambda row: ";".join(row.dropna().astype(str)), axis=1)
    # Remove all columns except Taxonomy
    tax = tax.drop(tax.columns.difference(["Taxonomy"]), axis=1)
    return tax

if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("-i", "--input", help="Input OTU table", required=True)
    args.add_argument("-f", "--fasta", help="Input FASTA file (OTUs)", required=True)
    args.add_argument("-t", "--taxonomy", help="Input taxonomy file (dadaist2 format)", required=True)
    args.add_argument("-o", "--output", help="Output OTU table", required=True)
    args.add_argument("-s", "--otutab-separator", help="Separator used in the OTU table", default='\t')
    args.add_argument("-k", "--index-name", help="Index name", default="#NAME")
    args.add_argument("--condense", help="Condense the taxonomy to a string", action="store_true")
    args.add_argument("--csv", help="Print CSV rather than TSV", action="store_true")
    args = args.parse_args()

    OUTPUT_SEP = "\t" if not args.csv else ","

    if not os.path.exists(args.input):
        print("Error: OTU table file not found: " + args.input)
        sys.exit(1)
    
    if not os.path.exists(args.taxonomy):
        print("Error: Taxonomy file not found: " + args.taxonomy)
        sys.exit(1)

    if not os.path.exists(args.fasta):
        print("Error: FASTA file not found: " + args.fasta)
        sys.exit(1)    
    # Read OTU table
    df = pd.read_csv(args.input, sep=args.otutab_separator, index_col=0)
    sampleNames = df.columns.values

    # Read taxonomy file
    taxonomy = loadTaxonomy(args.taxonomy, args.fasta)
    
    taxonomy_condensed = condenseTaxonomy(taxonomy)
    
    # Add the taxonomy joining df and taxonomy
    if args.condense:
        df = pd.concat([df, taxonomy_condensed], axis=1)
        taxonomy_condensed.index.name = "#TAXONOMY"
        taxonomy_condensed.to_csv(args.output + ".taxonomy.txt", sep=OUTPUT_SEP)
    else:
        df = pd.concat([df, taxonomy], axis=1)
        taxonomy.index.name = "#TAXONOMY"
        taxonomy.to_csv(args.output + ".taxonomy.txt", sep=OUTPUT_SEP)
    
    # Set index name
    df.index.name = args.index_name

    # Remove all rows with NaN values
    df = df.dropna(axis=0, how='any')
    
    df.to_csv(args.output, sep=OUTPUT_SEP)

    