#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

nextflow.enable.dsl = 2

/*
========================================================================================
                         MPRAflow
========================================================================================
MPRA Analysis Pipeline. Started 2019-07-29.
Library Association package

#### Homepage / Documentation
https://github.com/shendurelab/MPRAflow
#### Authors
Gracie Gordon <gracie.gordon@ucsf.edu>
Max Schubach <max.schubach@bihealth.de>
Sean Whalen <sean.whalen@gladstone.ucsf.edu>

### Modifications as of 2024-10-01
Adelaide Tovar <tovar@umich.edu>
Reworked code to be compatible with DSL 2
Created new option to not merge paired end reads and align separately
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     shendurelab/MPRAflow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run MPRA-nextflow -profile singularity,test

    Mandatory arguments:
      --fastq-insert                Full path to library association fastq for insert (must be surrounded with quotes)
      --fastq-bc                    Full path to library association fastq for bc (must be surrounded with quotes)
      --design-file                 Full path to fasta of ordered oligo sequences (must be surrounded with quotes)
      --design                      Base name of fasta
      --name                        Name of the association. Files will be named after this.

    Options:
      --fastq-insertPE              Full path to library association fastq for read2 if the library is paired end (must be surrounded with quotes)
      --merge-PE                     If providing read2, merge pairs with
      --min-cov                     minimum coverage of bc to count it (default 3)
      --min-frac                    minimum fraction of bc map to single insert (default 0.5)
      --mapq                        map quality (default 30)
      --baseq                       base quality (default 30)
      --cigar                       require exact match ex: 200M (default none)
      --outdir                      The output directory where the results will be saved and what will be used as a prefix (default outs)
      --split                       Number read entries per fastq chunk for faster processing (default: 2000000)
      --labels                      tsv with the oligo pool fasta and a group label (ex: positive_control) if no labels desired a file will be automatically generated

    Extras:
      --h, --help                   Print this help message
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.containsKey('h') || params.containsKey('help')){
    helpMessage()
    exit 0
}

//defaults
params.min_cov="2"
params.min_frac="0.5"
params.baseq="10"
params.mapq="1"
params.cigar="n"
params.split=2000000

// Validate inputs
if ( !params.containsKey("name") ){
    exit 1, "Please specify a name of this workflow using --name"
}

if ( params.containsKey("fastq-insert") ){
    params.fastq_insert_file = file(params['fastq-insert'])
    if( !params.fastq_insert_file.exists() ) exit 1, "Fastq insert file not found: ${params.fastq_insert_file}"
} else {
    exit 1, "Fastq insert file not specified with --fastq-insert"
}

if(params.containsKey("fastq-insertPE")){
    params.fastq_insertPE_file = file(params['fastq-insertPE'])
    if( !params.fastq_insertPE_file.exists() ) exit 1, "Fastq paired-end insert file not found: ${params.fastq_insertPE_file}"
} else {
  params.fastq_insertPE_file = null
}

if (params.containsKey("merge-PE")){
    params.mergePE = true
} else {
    params.mergePE = false
}

// Fastq barcode file in params.fastq_bc_file
if ( params.containsKey("fastq-bc")){
    params.fastq_bc_file = file(params['fastq-bc'])
    if( !params.fastq_bc_file.exists() ) exit 1, "Fastq barcode file not found: ${params.fastq_bc_file}"
} else {
    exit 1, "Fastq barcode file not specified with --fastq-bc"
}

// design file saved in params.design_file
if ( params.containsKey("design-file")){
    params.design_file = file(params['design-file'])
    if( !params.design_file.exists() ) exit 1, "Design file ${params.design_file} does not exist"
} else {
    exit 1, "Design file not specified with --design-file"
}

if ( params.containsKey("design")){
    params.design = file(params['design'])
} else {
    exit 1, "Design file name not given with --design"
}

// label file saved in params.label_file
if (params.containsKey("labels")){
    params.label_file = file(params['labels'])
    if (!params.label_file.exists()) exit 1, "Label file ${params.label_file} does not exist"
} else {
    params.label_file = null
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

MPRAflow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'MPRAflow'
summary['Pipeline Version'] = params.version
summary['Fastq insert']     = params.fastq_insert_file
summary['Fastq paired']     = params.fastq_insertPE_file
summary['Fastq barcode']    = params.fastq_bc_file
summary['Design fasta']     = params.design_file
summary['Design name']      = params.design
summary['Minimum BC cov']   = params.min_cov
summary['Map quality']      = params.mapq
summary['Base quality']     = params.baseq
summary['Cigar string']     = params.cigar
summary['Min % mapped']     = params.min_frac
summary['Output dir']       = params.outdir
summary['Run name']         = params.name
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Base directory']   = "$baseDir"
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile

//summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
//summary['Max CPUs']         = params.max_cpus
//summary['Max Time']         = params.max_time

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
* count fastq and bam length remove the illegal regex characters
* and make label file if missing
* contributions: Gracie Gordon & Max Schubach
*/

process 'count_bc' {

    tag 'count'
    publishDir "${params.outdir}/${params.name}", mode:'copy'

    input:
    path fastq_bc
    path design_file
    path labels

    output:
    path 'count_fastq.txt', emit: bc_output
    path "*.rmIllegalChars.fa", emit: fixed_design
    path "*.label.rmIllegalChars.txt", emit: fixed_label

    shell:
    def prefix = task.ext.prefix ?: "${design_file.baseName}"
    """
    ## Get rid of illegal regex characters
    awk '{gsub(/\\[/,"_")}1' $labels > t_new_label.txt
    awk '{gsub(/\\]/,"_")}1' t_new_label.txt > ${prefix}.label.rmIllegalChars.txt

    awk '{gsub(/\\[/,"_")}1' $design_file > t_new_design.txt
    awk '{gsub(/\\]/,"_")}1' t_new_design.txt > ${prefix}.rmIllegalChars.fa

    zcat $fastq_bc | wc -l  > count_fastq.txt
    """

}

/*
* count fastq and bam length remove the illegal regex characters
* and make design file
* contributions: Gracie Gordon & Max Schubach
*/

process 'count_bc_nolab' {

    tag 'count'
    publishDir "${params.outdir}/${params.name}", mode:'copy'

    input:
    path fastq_bc
    path design_file
    
    output:
    path 'count_fastq.txt', emit: bc_output
    path "*.rmIllegalChars.fa", emit: fixed_design
    path "*.label.rmIllegalChars.txt", emit: fixed_label

    shell:
    def prefix = task.ext.prefix ?: "${design_file.baseName}"
    """
    #CREATE LABEL FILE and remove illegal regex characters
    awk -F'\t' 'BEGIN {OFS = FS} NR%2==1 {print substr(\$1,2,length(\$1)),"na"}' $design_file > labels.txt
    awk '{gsub(/\\[/,"_")}1' labels.txt > t_new_label.txt
    awk '{gsub(/\\]/,"_")}1' t_new_label.txt > ${prefix}.label.rmIllegalChars.txt

    awk '{gsub(/\\[/,"_")}1' $design_file | \\
    awk '{gsub(/\\]/,"_")}1' | \\
    sed 's/\\r//g' > ${prefix}.rmIllegalChars.fa

    zcat $fastq_bc | wc -l  > count_fastq.txt
    """

}

/*
* STEP 1: Align
* Process 1A: create BWA reference
* contributions: Gracie Gordon
*/

process 'create_BWA_ref' {

    tag "make ref"
    conda 'conf/mpraflow_py36.yml'
    cpus 4
    memory '64 GB'
    time '4h'
    queue 'standard'

    input:
    path(fixed_design)

    output:
    path("${fixed_design}.*"), emit: reference

    shell:
    """
    bwa index -a is ${fixed_design}
    samtools faidx ${fixed_design}
    java -jar $PICARDLIB/picard.jar CreateSequenceDictionary --REFERENCE ${fixed_design} --OUTPUT ${fixed_design}.dict
    """

}

/*
*Process 1B: merge Paired end reads
* contributions: Gracie Gordon
*/
process 'PE_merge' {

    tag 'merge'
    conda 'conf/mpraflow_py36.yml'

    input:
    tuple(val(chunk_id), path(R1_chunk), path(R3_chunk))
    val(params.name)

    output:
    path("${params.name}.${chunk_id}.merged.*.fastq"), emit: mergedPE

    shell:
    """
    fastq-join $R1_chunk $R3_chunk -o ${params.name}.${chunk_id}.merged.%.fastq
    """

}

/*
* Process 1C: align with BWA
* contributions: Gracie Gordon
*/

//paired ends
process 'align_BWA_PE' {

    tag "align"
    cpus 10
    memory '16 GB'
    time '12h'
    queue 'standard'

    conda 'conf/mpraflow_py36.yml'

    input:
    path(fixed_design)
    path(mergedPE)
    val(params.name)
    path(reference)

    output:
    path('*.sorted.bam'), emit: s_bam
    path('*.count_bam.txt'), emit: bam_ch

    shell:
    """
    bwa mem ${fixed_design} ${params.name}.${chunk_id}.merged.join.fastq | samtools sort -o ${params.name}.${chunk_id}.sorted.bam

    echo 'bam made'

    samtools view ${params.name}.${chunk_id}.sorted.bam | wc -l > ${params.name}.${chunk_id}.count_bam.txt
    """

}

//paired end, not merged
process 'align_BWA_PE_separate' {

    tag "align"
    cpus 10
    memory '16 GB'
    time '12h'
    queue 'standard'

    conda 'conf/mpraflow_py36.yml'

    input:
    path(fixed_design)
    tuple val(chunk_id), path(R1_chunk), path(R3_chunk)
    val(params.name)
    path(reference)

    output:
    path('*.sorted.bam'), emit: s_bam
    path('*.count_bam.txt'), emit: bam_ch

    shell:
    """
    bwa mem ${fixed_design} $R1_chunk $R3_chunk | samtools sort -o ${params.name}.${chunk_id}.sorted.bam

    echo 'bam made'

    samtools view ${params.name}.${chunk_id}.sorted.bam | wc -l > ${params.name}.${chunk_id}.count_bam.txt
    """

}

//single end
process 'align_BWA_S' {

    tag "align"
    cpus 10
    memory '16 GB'
    time '12h'
    queue 'standard'

    conda 'conf/mpraflow_py36.yml'

    input:
    path(fixed_design)
    tuple val(chunk_id), path(R1_ch)
    val(params.name)
    path(reference)

    output:
    path('*.sorted.bam'), emit: s_bam
    path('*.count_bam.txt')

    shell:
    """
    bwa mem ${fixed_design} $R1_ch | samtools sort - o ${params.name}.${chunk_id}.sorted.bam

    echo 'bam made'

    samtools view ${params.name}.${chunk_id}.sorted.bam | wc -l > ${params.name}.${chunk_id}.count_bam.txt
    """

}

/*
*COLLCT FASTQ CHUNCKS
*/

process 'collect_chunks'{

    cpus 2
    memory '160 GB'
    time '24h'
    queue 'standard'

    conda 'conf/mpraflow_py36.yml'

    input:
    path s_bam

    output:
    path('s_merged.bam'), emit: s_merge
    path('count_merged.txt'), emit: ch_merge

    shell:
    """
    #collect sorted bams into one file
    samtools merge -f all.bam ${s_bam.join(' ')}
    samtools sort all.bam -o s_merged.bam

    #collect bam counts into one file
    samtools view s_merged.bam | wc -l > count_merged.txt
    """

}

/*
* Assign barcodes to element sequences
* contributions: Sean Whalen
*/

process 'map_element_barcodes' {

    tag "assign"
    cpus 10
    memory '240 GB'
    time '24h'
    queue 'largemem'
    publishDir "${params.outdir}/${params.name}", mode:'copy'

    conda 'conf/mpraflow_py36.yml'

    input:
    val(params.name)
    val(params.mapq)
    val(params.baseq)
    val(params.cigar)
    path(fastq_bc)
    path(bc_output)
    path(s_merge)
    path(ch_merge)

    output:
    path("${params.name}_coords_to_barcodes.pickle"), emit: map_ch
    path("${params.name}_barcodes_per_candidate-no_repeats-no_jackpots.feather"), emit: count_table_ch
    path("${params.name}_barcode_counts.pickle"), emit: pickled

    shell:
    """
    echo "test assign inputs"
    echo ${params.mapq}
    echo ${params.baseq}
    echo $fastq_bc
    echo ${bc_output}
    echo ${ch_merge}

    python ${"$baseDir"}/src/nf_ori_map_barcodes.py ${"$baseDir"} ${fastq_bc} ${bc_output} ${ch_merge} ${s_merge} ${params.name} ${params.mapq} ${params.baseq} ${params.cigar}
    """

}

/*
* Filter barcodes for minimum coverage and unique mapping
* contributions: Gracie Gordon
*/

process 'filter_barcodes' {

    tag "filter"
    cpus 10
    memory '80 GB'
    time '24h'
    queue 'standard'
    publishDir "${params.outdir}/${params.name}", mode:'copy'

    conda 'conf/mpraflow_py36.yml'

    input:
    val(params.min_cov)
    val(params.min_frac)
    val(params.name)
    path(map_ch)
    path(count_table_ch)
    path(fixed_label)

    output:
    path("${params.name}_filtered_coords_to_barcodes.pickle")
    path("${params.name}_original_counts.png")
    path("original_count_summary.txt")
    path("${params.name}_filtered_counts.png")
    path("filtered_count_summary.txt")

    shell:
    """
    python ${"$baseDir"}/src/nf_filter_barcodes.py ${params.name} ${map_ch} ${count_table_ch} ${params.min_cov} ${params.min_frac} ${fixed_label}
    """

}

/*
Define workflow
*/

workflow {

    R1_ch = Channel.fromPath(params.fastq_insert_file)
                    .splitFastq(by: params.split, decompress: true, file: true)
                    .map { file ->
                        def chunk_id = file.name.tokenize('.')[2]
                        tuple(chunk_id, file)    
                    }

    if (params.fastq_insertPE_file != null) {
        R3_ch = Channel.fromPath(params.fastq_insertPE_file)
                    .splitFastq(by: params.split, decompress: true, file: true)
                    .map { file ->
                        def chunk_id = file.name.tokenize('.')[2]
                        tuple(chunk_id, file)
                    }

        paired_ch = R1_ch.join(R3_ch)
                            .map { chunk_id, file1, file2 ->
                                tuple(chunk_id, file1, file2)
                            }

    }

    bc_ch = Channel.value(params.fastq_bc_file)
    design_ch = Channel.value(params.design_file)

    if (params.label_file) {
        label_ch = Channel.value(params.label_file)
        count_bc(bc_ch, design_ch, label_ch)
        bc_out_ch = count_bc.out.bc_output.collect()
        fixed_design_ch = count_bc.out.fixed_design.collect()
        fixed_label_ch = count_bc.out.fixed_label.collect()
    } else {
        count_bc_nolab(bc_ch, design_ch)
        bc_out_ch = count_bc_nolab.out.bc_output.collect()
        fixed_design_ch = count_bc_nolab.out.fixed_design.collect()
        fixed_label_ch = count_bc_nolab.out.fixed_label.collect()
    }

    create_BWA_ref(fixed_design_ch)
    ref_ch = create_BWA_ref.out.reference.collect()

    if (params.fastq_insertPE_file != null && params.mergePE) {
        merged_ch = PE_merge(paired_ch, params.name)
        (sorted_bam, counts) = align_BWA_PE(fixed_design_ch, merged_ch, params.name, ref_ch)
    } else if (params.fastq_insertPE_file != null) {
        (sorted_bam, counts) = align_BWA_PE_separate(fixed_design_ch, paired_ch, params.name, ref_ch)
    } else {
        (sorted_bam, counts) = align_BWA_S(fixed_design_ch, R1_ch, params.name, ref_ch)
    }

    (full_bam, full_count) = collect_chunks(sorted_bam.collect())
    (map_ch, count_table_ch, pickled) = map_element_barcodes(params.name, params.mapq, params.baseq, params.cigar, bc_ch, bc_out_ch, full_count, full_bam)
    filter_barcodes(params.min_cov, params.min_frac, params.name, map_ch, count_table_ch, fixed_label_ch)

}
