#!/usr/bin/env python3
"""
#NAME,Files
A01,"A01_S0_L001_R1_001.fastq.gz,A01_S0_L001_R2_001.fastq.gz"
A02,"A02_S0_L001_R1_001.fastq.gz,A02_S0_L001_R2_001.fastq.gz"
F99,"F99_S0_L001_R1_001.fastq.gz,F99_S0_L001_R2_001.fastq.gz"


#NAME,A01,A02,F99
ASV1,1263,1544,1341
ASV2,116,89,57
ASV3,100,21,19
ASV4,77,13,11
ASV5,69,18,12
ASV6,29,40,0

#TAXONOMY,Domain,Phylum,Class,Order,Family,Genus,Species
ASV1,d__Bacteria,p__Proteobacteria,c__Betaproteobacteria,o__Burkholderiales,f__Burkholderiaceae,g__Ralstonia,s__
ASV2,d__Bacteria,p__,c__,o__,f__,g__,s__
ASV3,d__Bacteria,p__Firmicutes,c__Bacilli,o__Lactobacillales,f__Streptococcaceae,g__Lactococcus,s__
ASV4,d__Bacteria,p__Firmicutes,c__Bacilli,o__Lactobacillales,f__Streptococcaceae,g__Lactococcus,s__
ASV5,d__Bacteria,p__Firmicutes,c__Bacilli,o__Lactobacillales,f__Carnobacteriaceae,g__Carnobacterium,s__
"""
import shutil
import os
import sys
import argparse
from random import randint
__VERSION__ = "0.1"


class Features:
    def __init__(self, fasta):
        self.fasta_file = fasta
        self.features = []
        self.sequences = {}
        self.comments = {}
        self.readFeatures()

    def readFeatures(self):
        for name, comm, seq in self._read_fasta(self.fasta_file):
            self.sequences[name] = seq
            self.comments[name] = comm
            self.features.append(name)

    def _read_fasta(self, path):
        if path.endswith(".gz"):
            import gzip
        
        seqName = None
        seqComment = None
        with (gzip.open if path.endswith('.gz') else open)(path, 'rt') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    if seqName is not None:
                        yield seqName, seqComment, sequence
                    seqName = line[1:].split()[0]
                    seqComment = line[1:].split()[1:] if len(line[1:].split()) > 1 else ""
                    sequence = ""
                    
                else:
                    sequence += line.strip()
        yield seqName, seqComment, sequence

class Metadata:
    def __init__(self, metadata_file, separator, header_char, id_field, quote_char="\""):
        self.metadata_file = metadata_file
        self.metadata = {}
        self.separator = separator
        self.header_char = header_char
        self.id_field = id_field
        self.id_field_index = 0
        self.fields = []
        self.samples = []
        self.quote = quote_char
        self.load()
    
    def load(self):
        if self.metadata_file is None:
            return
        with open(self.metadata_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(self.header_char):
                    fields = line[len(self.header_char):].split(self.separator)
                    if self.id_field not in fields:
                        raise Exception(f"ID field not found in metadata file: {self.id_field} not in {fields}")
                    else:
                        # Determine the column of ID fields
                        self.id_field_index = fields.index(self.id_field)
                        self.fields = fields
                    continue
                if line == "":
                    continue
                fields = line.split(self.separator)
                id = fields[self.id_field_index]
                self.metadata[id] = fields
                self.samples.append(id)
    
    def csv(self):
        # 
        other_fields = [f for f in self.fields if f != self.id_field]
        # If a field contains ",", surround it by quotes

                
        csv = "#NAME," + ",".join(other_fields) + "\n"
        for sample in self.metadata:
            fields = self.metadata[sample]
            fields = [sample] + fields[self.id_field_index + 1:]
            for i in range(len(fields)):
                if fields[i].find(",") != -1: 
                    fields[i] = f'{self.quote}{fields[i]}{self.quote}'

            csv += ",".join(fields) + "\n"

        return csv

class FeatureTable:
    def __init__(self, feature_table_file, separator, header_char, id_field, quote_char="\""):
        self.feature_table_file = feature_table_file
        self.feature_table = {}
        self.separator = separator
        self.header_char = header_char
        self.id_field = id_field
        self.id_field_index = 0
        self.samples = []
        self.features = []
        self.quote = quote_char
        self.load()
    
    def load(self):
       
        with open(self.feature_table_file, "r") as f:
            for line in f:
                line = line.strip()
                # Parse header
                if line.startswith(self.header_char):
                    fields = line[len(self.header_char):].split(self.separator)
                    if self.id_field not in fields:
                        raise Exception(f"ID field not found in feature table file: {self.id_field} not in {fields}")
                    else:
                        # Determine the column of ID fields
                        self.id_field_index = fields.index(self.id_field)
                        self.samples = fields[1:]
                    continue
                if line == "":
                    continue
                fields = line.split(self.separator)
                id = fields[self.id_field_index]
                
                self.feature_table[id] = fields[1:]
                self.features.append(id)
    
    def csv(self):
        # 
        other_fields = [f for f in self.samples if f != self.id_field]
        # If a field contains ",", surround it by quotes
        csv = "#NAME," + ",".join(other_fields) + "\n"
        for sample in self.feature_table:
            fields = self.feature_table[sample]
            fields = [sample] + fields
            for i in range(len(fields)):
                if fields[i].find(",") != -1: 
                    fields[i] = f'{self.quote}{fields[i]}{self.quote}'

            csv += ",".join(fields) + "\n"

        return csv

def MakeMetadata(FeatureTable):
    metadata = Metadata(None, "\t", "#", "SampleID")
    metadata.samples = FeatureTable.samples
    metadata.fields = ["RandomGroup"]
    for sample in FeatureTable.samples:
        Letter = "A" if randint(0,1) == 0 else "B"
        metadata.metadata[sample] = [sample, Letter]
    return metadata

def loadDecipherTaxonomy(path, Sequences):
    """
    Kingdom Phylum Class Order Family Genus
    1 Bacteria Firmicutes Negativicutes Veillonellales-Selenomonadales Veillonellaceae Veillonella
    2 Bacteria Firmicutes Negativicutes Veillonellales-Selenomonadales Veillonellaceae Veillonella
    3 Bacteria Firmicutes Negativicutes Veillonellales-Selenomonadales Veillonellaceae Veillonella
    4 Bacteria Firmicutes Negativicutes Veillonellales-Selenomonadales Veillonellaceae Veillonella
    5 Bacteria Bacteroidota Bacteroidia Bacteroidales Prevotellaceae Prevotella
    6 Bacteria Firmicutes Negativicutes Veillonellales-Selenomonadales Veillonellaceae Veillonella
    7 Bacteria Actinobacteriota Actinobacteria Micrococcales Brevibacteriaceae Brevibacterium
    8 Bacteria Halanaerobiaeota Halanaerobiia Halanaerobiales Halanaerobiaceae Halanaerobium
    9 Bacteria Firmicutes Negativicutes Veillonellales-Selenomonadales Veillonellaceae Veillonella
    """
    taxonomy = {}
    with open(path, "r") as f:
        c = 0
        
        for line in f:
            line = line.strip().split()
            
            if c >= 1:
                seq = Sequences.features[c-1]
                taxonomy[seq] = line[1:]

            c += 1
    return taxonomy

def taxonomyToCsv(taxonomy):
    ranks = "Domain,Phylum,Class,Order,Family,Genus,Species".split(",")
    csv = "#NAME," + ",".join(ranks) + "\n"
    for feature in taxonomy:
        fields =  taxonomy[feature]
        for i in range(len(ranks)):
            # add "" to fields if shorter than ranks
            if len(fields) < i+1:
                fields.append("")
            if fields[i].find(",") != -1: 
                fields[i] = f'"{fields[i]}"'
            fields[i] = "" if fields[i] == "NA" else fields[i]
            fields[i] = ranks[i][:1].lower() + "__" + fields[i]

        csv += feature + "," + ",".join(fields) + "\n"

    return csv
if __name__ == "__main__":
    args = argparse.ArgumentParser("Export UFLOW to MicrobiomeAnalyst and Phyloseq")
    args.add_argument("-i", "--feature-table", help="Input OTU table", required=True)
    args.add_argument("-f", "--fasta", help="Input FASTA file (OTUs)", required=True)
    args.add_argument("-t", "--taxonomy", help="Input taxonomy file (dadaist2 format)", required=True)
    args.add_argument("-m", "--metadata", help="Metadata file", required=False)
    args.add_argument("--tree", help="Tree file (copied)", required=False)
    args.add_argument("-o", "--output", help="Output directory [default: %(default)s]", default="MicrobiomeAnalyst")
    args.add_argument("--verbose", help="Verbose output", action="store_true")

    # Add section for parsers
    parserargs = args.add_argument_group("Parsing options")
    parserargs.add_argument("--meta-sep", help="Metadata separator [default: tab]", default="\t")
    parserargs.add_argument("--meta-header", help="Metadata header char [default: %(default)s]", default="#")
    parserargs.add_argument("--meta-id", help="Metadata header identifier [default: %(default)s]", default="SampleID")
    parserargs.add_argument("--table-sep", help="Feature table separator [default: tab]", default="\t")
    parserargs.add_argument("--table-header", help="Feature table header char [default: %(default)s]", default="#")
    parserargs.add_argument("--table-id", help="Feature table header identifier [default: %(default)s]", default="OTU ID")
    
    args.add_argument("--version", action="version", version="%(prog)s " + __VERSION__)
    args = args.parse_args()

    if not os.path.exists(args.output):
        if args.verbose:
            print("Creating output directory: {}".format(args.output))
        os.makedirs(args.output)
    
    # Check input files existance
    if not os.path.exists(args.feature_table):
        print("ERROR: Feature table file does not exist: {}".format(args.feature_table))
        sys.exit(1)
    if not os.path.exists(args.fasta):
        print("ERROR: FASTA file does not exist: {}".format(args.fasta))
        sys.exit(1)
    if not os.path.exists(args.taxonomy):
        print("ERROR: Taxonomy file does not exist: {}".format(args.taxonomy))
        sys.exit(1)
    
    if args.metadata is not None and not os.path.exists(args.metadata):
        print("ERROR: Metadata file does not exist: {}".format(args.metadata))
        sys.exit(1)


    # Load feature table
    FeatTable = FeatureTable(args.feature_table, args.table_sep, args.table_header, args.table_id)

    # Load or make metadata
    metadata = None
    if args.metadata is not None:
        metadata = Metadata(args.metadata, args.meta_sep, args.meta_header, args.meta_id)
    else:        
        metadata = MakeMetadata(FeatTable)
  
    RepSeqs = Features(args.fasta)

    # Some checks?
    if len(FeatTable.samples) != len(metadata.samples):
        print("ERROR: Number of samples in feature table and metadata do not match")
        sys.exit(1)

    if len(FeatTable.features) != len(RepSeqs.features):
        print("ERROR: Number of features in feature table and FASTA do not match")
        sys.exit(1)

    # Compare FeatTable.features != RepSeqs.features:
    for f in FeatTable.features:
        if f not in RepSeqs.features:
            print("ERROR: Feature {} not found in FASTA".format(f))
            sys.exit(1)
    tax = loadDecipherTaxonomy(args.taxonomy, RepSeqs)

    # Make output files
    if args.tree is not None:
        dest = os.path.join(args.output, "rep-seqs.tree")
        shutil.copy(args.tree, dest)
    with open(os.path.join(args.output, "metadata.csv"), "w") as f:
        f.write(metadata.csv())

    with open(os.path.join(args.output, "table.csv"), "w") as f:
        f.write(FeatTable.csv())
    with open(os.path.join(args.output, "taxonomy.csv"), "w") as f:
        print(taxonomyToCsv(tax), file=f)
