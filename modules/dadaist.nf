process TAX {
    tag "${task.cpus}"
    label 'process_medium'
    publishDir "$params.outdir/", 
        mode: 'copy'

    input:
    path("seqs.fa")
    path("DB")

    
    output:
    path("taxonomy"), optional: true
    
    script:
    """
    parseDb.py -i DB
    dadaist2-assigntax -i seqs.fa --outdir taxonomy -t ${task.cpus} --reference DB.*
    """    
}