/*  Input parameters   */
nextflow.enable.dsl = 2
params.dir        = "$baseDir/reads/"
params.pattern    = "_R{1,2}.fastq.gz"
def reads         = params.dir + "/*" + params.pattern
def db            = params.dbdir + "/" + params.db

params.outdir     = "$baseDir/uflow"
params.its        = false
params.its_region = "ITS1"
params.forward    = "CTTGGTCATTTAGAGGAAGTAA"
params.reverse    = "GCTGCGTTCTTCATCGATGC"  

params.skip_uncross = false
      
// prints to the screen and to the log
log.info """
         GMH Amplicon Analysis (version 1.3)
         ===================================
         input reads  : ${reads}
         database     : ${db}
         outdir       : ${params.outdir}
         primers      : ${params.forward}:${params.reverse}
         specs        : ${params.max_cpus} cores, ${params.max_memory} GB memory
         """
         .stripIndent()

/* 
   check reference path exists 
*/

def dbPath = file(db, checkIfExists: true)
 /*    Modules  */
include { CUTADAPT; ITSX; MERGE; RELABEL; FASTP; FILT; DEREP; UNOISE;  } from './modules/amplicon'
include {  OTUTABLE; JOINTAB; NORM; ADDTAX; UNCROSS; TABSTATS; OCTAVE; ALPHA; BETA; TRIM } from './modules/otutab'
include { TAX } from './modules/dadaist'
reads = Channel
        .fromFilePairs(reads, checkIfExists: true)

workflow CLEAN {
  take:
    reads
    fwd_primer
    rev_primer

  main:
    // Takes FASTQ Pe, returns FASTQ Pe
    FASTP(reads)
    RELABEL(FASTP.out.reads)
    CUTADAPT(RELABEL.out, fwd_primer, rev_primer )

  emit:
    CUTADAPT.out   
}


workflow {
  // Discard samples not passing the min reads filter
  CLEAN(reads, params.forward, params.reverse)


  if (params.its == true) {
      READS = ITSX(CLEAN.out, params.its_region)
  } else {
      READS = MERGE(CLEAN.out)
  }
  FILT(READS)
  DEREP(FILT.out.map{it -> it[1]}.collect())
  UNOISE(DEREP.out)
  TAX(UNOISE.out, dbPath)
  OTUTABLE(READS, UNOISE.out)

  JOINTAB(OTUTABLE.out.map{it -> it[1]}.collect())
  UNCROSS(JOINTAB.out, params.skip_uncross)
  
  TRIM(UNCROSS.out.table)

  TABSTATS(TRIM.out)
  ALPHA(TRIM.out)
  BETA(TRIM.out)
  
  OCTAVE(TRIM.out, UNOISE.out)
  NORM(TRIM.out)
  ADDTAX(NORM.out, TAX.out, UNOISE.out)
  
  //TRACKFILES(FASTP.out.json.mix( KRAKEN2_HOST.out.txt, CONTAMLOG ).collect() )
  //MULTIQC( FASTP.out.json.mix( KRAKEN2_REPORT.out, TRACKFILES.out ).collect() )
}