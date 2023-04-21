NAME=uflow
mamba create -n $NAME -y -c conda-forge -c bioconda \
  nextflow pigz "seqfu>=1.9" "multiqc>1.9" "fastp" \
  "vsearch=2.18" "clustalo" "fasttree" "itsxpress" \
  "pandas>1.0" "cutadapt" "dadaist2"

