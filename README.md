# uflow

USEARCH based workflow. See [dadaist2](https://github.com/quadram-institute-bioscience/dadaist2) for an Open Source workflow.

## Instructions

1. Generate a singularity container using the definition from [`deps`](./deps/)
2. Download Dadaist2 databases and use that directory as `--dbdir`
3. Run the workflow selecting or creating a profile:
```bash
nextflow run uflow/ --input DIR --outdir OUT -profile bact
```

Parameters:
* `--outdir` output directory
* `--dir` path to the input directory
* `--pattern` expansion to get the read pairs (default = "_R{1,2}.fastq.gz")
* `--dbdir` path to the database directory (hint: put it in the profile)
* `--db` actual database (file in the dbdir) to be used
* `--its` enable ITS processing
* `--its_region` ITS region (default = "ITS1")
* `--forward` sequence of the forward primer (default = "CTTGGTCATTTAGAGGAAGTAA")
* `--reverse`  sequence of the reverse primer (default = "GCTGCGTTCTTCATCGATGC")
