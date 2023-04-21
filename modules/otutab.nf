process JOINTAB {
    label 'process_low'
    publishDir "$params.outdir/", 
        mode: 'copy'
    input:
    path("*.tab")
    
    output:
    path("rawtable.tsv")

    script:
    """
    mergeTables.py -v --sort --sample-from-header --otuid "#OTU ID" -o rawtable.tsv *.tab
    """        
}

process NORM {

    label 'process_low'

    input:
    path("raw.tab")
    
    output:
    path("freq_table.tsv")

    script:
    """
    normalizeOtutable.py  -i raw.tab  -o freq_table.tsv 
    """            
}

process OTUTABLE {
    tag "filter $sample_id"
    label 'parallel'
    
    input:
    tuple val(sample_id), path(reads) 
    path("asv.fasta")
    
    
    output:
    tuple val(sample_id), path("${sample_id}.tab")
    
    script:
    """
    usearch -otutab ${reads[0]} -otutabout table.tmp -otus asv.fasta -threads ${task.cpus}
    # Fix sample name: remplace the second column of the fist line by the sample name
    sed  "1s/.*/#OTU\t${sample_id}/" table.tmp >  ${sample_id}.tab
    """    
}
process TABSTATS {
    label 'process_medium'
    publishDir "$params.outdir/info/",
        mode: 'copy'    

    input:
    path("otutab.txt")
    
    
    output:
    path("otutable-stats.txt"), emit: html, optional: true

    """
    usearch -otutab_stats otutab.txt -output otutable-stats.txt
    """
}

process TRIM {
    label 'process_medium' 

    input:
    path("otutab.txt")
    
    
    output:
    path("otutable.tsv")

    """
    usearch -otutab_trim otutab.txt \
      -min_sample_size 1000 -min_freq 0.0005 -min_otu_size 10 \
      -output otutable.tsv
    """
}

process OCTAVE {
    label 'process_high'
    publishDir "$params.outdir/info/",
        mode: 'copy'    

    input:
    path("otutab.txt")
    path("seqs.fa")
    
    
    output:
    path("octave.*"),  optional: true

    """
    set +e
    usearch -calc_distmx seqs.fa -tabbedout distmx.txt -maxdist 0.2 -termdist 0.3 -threads ${task.cpus} || touch fail.txt
    if [[ ! -e fail.txt ]]; then
        usearch -otutab_octave otutab.txt -distmxin distmx.txt \
            -htmlout octave.html -svgout octave.svg || true
    fi
    """
} 
process ALPHA {
    label 'process_medium'
    publishDir "$params.outdir/diversity/",
        mode: 'copy'    

    input:
    path("otutab.txt")

    
    output:
    path("alpha.txt"),  optional: true

    """
    usearch -alpha_div otutab.txt -output alpha.txt
    """
} 

process BETA {
    label 'process_medium'
    publishDir "$params.outdir/diversity/",
        mode: 'copy'    

    input:
    path("otutab.tab")

    
    output:
    path("*.{txt,tree}"),  optional: true

    """
    usearch -beta_div otutab.tab -metrics jaccard,bray_curtis
    """
} 

process UNCROSS {
    tag "$skip"
    label 'process_medium'
    publishDir "$params.outdir/info/",
        pattern: '*.html',
        mode: 'copy'    

    input:
    path("otutab.txt")
    val(skip)
    
    output:
    path("otutab.tsv"), emit: table
    path("crosstalk.html"), emit: html, optional: true

    script:
    if (skip == false)
        """
        set +e
        usearch -otutab_xtalk otutab.txt -report xtalk_report.txt \
            -htmlout crosstalk.html -otutabout otutab.tsv || cp otutab.txt otutab.tsv
        """
   else
       """
       cp  otutab.txt otutab.tsv
       """
}

process ADDTAX {
    tag "filter $sample_id"
    label 'parallel'
    publishDir "$params.outdir/", 
        mode: 'copy'
    
    input:
    path("table.tsv")
    path("taxonomy")
    path("asv.fa")
    
    
    output:
    path("table_relative_tax.tsv")
    
    script:
    """
    addTaxonomy.py -t taxonomy/taxonomy.tsv -f asv.fa -i table.tsv -o table_relative_tax.tsv
    """    
}


