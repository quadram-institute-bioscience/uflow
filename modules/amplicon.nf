
process CUTADAPT {
    tag "primers $sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads) 
    val(fwd)
    val(rev)
    
    output:
    tuple val(sample_id), path("filtered/${sample_id}_R*.fastq.gz")
    
    script:
    """
    FORx=\$(seqfu rc ${fwd})
    REVx=\$(seqfu rc ${rev})
    mkdir -p filtered
    cutadapt -a ${fwd}...\${REVx} -A ${rev}...\${FORx}  --discard-untrimmed  -j ${task.cpus} \
        -o filtered/${sample_id}_R1.fastq.gz -p filtered/${sample_id}_R2.fastq.gz ${reads[0]} ${reads[1]}
    """
}


process ITSX {
    tag "ITSxpress $sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads) 
    val(region)
    
    output:
    tuple val(sample_id), path("${sample_id}.fastq")
    
    script:
    """
    itsxpress --region ${region} --threads ${task.cpus} \
     --fastq ${reads[0]} --fastq2 ${reads[1]} --outfile ${sample_id}.fastq
    
    """
}
process FASTP {
    tag "filter $sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads) 

    
    output:
    tuple val(sample_id), path("filt/${sample_id}_R*.fastq.gz"), emit: reads
    tuple val(sample_id), path("filt/${sample_id}_R*.fastq.gz"), emit: json
    script:
    """
    mkdir -p filt
    fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} -o filt/${sample_id}_R1.fastq.gz -O filt/${sample_id}_R2.fastq.gz \
      --detect_adapter_for_pe --n_base_limit 1 --length_required 50 -j ${sample_id}.fastp.json
    """
}

process RELABEL {
    tag "rename $sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(reads) 

    
    output:
    tuple val(sample_id), path("out/${sample_id}_R*.fastq.gz")

    script:
    """
    mkdir -p out
    seqfu cat -p ${sample_id}. -z ${reads[0]}  | gzip -c > out/${sample_id}_R1.fastq.gz
    seqfu cat -p ${sample_id}. -z ${reads[1]}  | gzip -c > out/${sample_id}_R2.fastq.gz
    """
}

process MERGE {
    tag "filter $sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads) 

    
    output:
    tuple val(sample_id), path("${sample_id}.fastq"), emit: reads optional true
    
    script:
    """
    #usearch -fastq_mergepairs ${reads[0]} -reverse ${reads[1]} -fastqout ${sample_id}.fastq -threads ${task.cpus} \
    #  -fastq_pctid 85 -fastq_maxdiffs 12
    
    vsearch --threads ${task.cpus} --fastq_mergepairs ${reads[0]} --reverse ${reads[1]} --fastqout ${sample_id}.fastq \
      --fastq_maxdiffs 14 --fastq_maxdiffpct 85.0 --fastq_maxns 0 --fastq_minmergelen 303

    # if empty remove file
    if [ ! -s "${sample_id}.fastq" ]; then
        echo "removing empty ${sample_id}.fastq"
        rm "${sample_id}.fastq"
    fi
    """
}

process FILT {
    tag "filter $sample_id"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads) 

    
    output:
    tuple val(sample_id), path("${sample_id}.fa")
    
    script:
    """
    usearch -fastq_filter ${reads[0]} -fastq_maxee 1.5 -relabel filt. -fastaout ${sample_id}.fa
    """    
}


process DEREP {
    tag "filter $sample_id"
    label 'process_medium'

    input:
    path("*")

    
    output:
    path("uniques.fasta")
    
    script:
    """
    seqfu cat *.fa | seqfu derep > uniques.fasta
    """    
}
process UNOISE {
    tag "filter $sample_id"
    label 'process_medium'
    publishDir "$params.outdir/", 
        mode: 'copy'

    input:
    path("uniques.fasta")

    
    output:
    path("asv.fasta")
    
    script:
    """
    usearch -unoise3 uniques.fasta -zotus asv.fasta
    """    
}
