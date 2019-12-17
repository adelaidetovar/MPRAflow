#!/usr/bin/env nextflow

params.version="2.0"
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
      --fastq_insert                Full path to library association fastq for insert (must be surrounded with quotes)
      --fastq_bc                    Full path to library association fastq for bc (must be surrounded with quotes)
      --design                      Full path to fasta of ordered oligo sequences (must be surrounded with quotes) 

    Options:
      --fastq_insertPE              Full path to library association fastq for read2 if the library is paired end (must be surrounded with quotes)
      --min_cov                     minimum coverage of bc to count it (default 3)
      --min_frac                    minimum fraction of bc map to single insert (default 0.5)
      --mapq                        map quality (default 30)
      --baseq                       base quality (default 30)
      --cigar                       require exact match ex: 200M (default none) 
      --outdir                      The output directory where the results will be saved and what will be used as a prefix (default outs)
      -w                            specific name for work directory (default: work)  
      -with-timeline                Create html file showing processing times
      --split                       number read entries per fastq chunk for faster processing (default: 2000000)  
      --labels                      tsv with the oligo pool fasta and a group label (ex: positive_control) if no labels desired a file will be automatically generated  

    Extras:
      --email                        Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false
output_docs = file("$baseDir/docs/output.md")

//defaults
params.min_cov="3"
params.min_frac="0.5"
params.baseq="30"
params.mapq="30"
params.cigar="n"
params.outdir="outs"
params.out=params.outdir
params.nf_required_version="19.07.0"
params.fastq_insertPE=0
params.split=2000000
params.labels=0


// Validate inputs
if ( params.fastq_insert ){
    fastq_insert = file(params.fastq_insert)
    if( !fastq_insert.exists() ) exit 1, "Fastq insert file not found: ${params.fastq_insert}"
}

if(params.fastq_insertPE != 0){
    fastq_insertPE = file(params.fastq_insertPE)
}

if ( params.fastq_bc ){
    fastq_bc = file(params.fastq_bc)
    if( !fastq_bc.exists() ) exit 1, "Fastq barcode file not found: ${params.fastq_bc}"
}

if ( params.design ){
    design = file(params.design)
    if( !design.exists() ) exit 1, "Fasta oligo design file not found: ${params.design}"
}

if (params.labels !=0){
    labels=file(params.labels)
    if (!labels.exists()) exit 1, "label file not specified ${params.labels}"
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
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Fastq insert']     = params.fastq_insert
summary['fastq paired']     = params.fastq_insertPE
summary['Fastq barcode']    = params.fastq_bc
summary['design fasta']     = params.design
summary['minimum BC cov']   = params.min_cov
summary['map quality']      = params.mapq
summary['base quality']     = params.baseq
summary['cigar string']     = params.cigar
summary['min % mapped']     = params.min_frac
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['base directory']   = "$baseDir"
summary['Working dir']      = workflow.workDir
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile

//summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
//summary['Max CPUs']         = params.max_cpus
//summary['Max Time']         = params.max_time


if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}




/*
* count fastq and bam length remove the illegal regex characters
* and make label file if missing
*/

if (params.labels != 0){
    process 'count_bc' {
        tag 'count'
        label 'shorttime'
        publishDir params.outdir, mode:'copy'
        
        input:
        file(fastq_bc) from fastq_bc
        file(design) from design
        //file(params.condaloc)
        file(label) from labels
        
        output:
        file 'count_fastq.txt' into bc_ch
        file "new_label.txt" into fixed_label
        file "new_design.fa" into fixed_design
        """
        #!/bin/bash
        cv=\$(which conda)
        cv1=\$(dirname "\$cv")
        cv2=\$(dirname "\$cv1")
        cv3=\${cv2}"/bin/activate"
        echo \$cv3
        source \$cv3 mpraflow_py36
        
        ## Get rid of illegal regex characters   
        awk '{gsub(/\\[/,"_")}1' $label > t_new_label.txt 
        awk '{gsub(/\\]/,"_")}1' t_new_label.txt > new_label.txt
        
        awk '{gsub(/\\[/,"_")}1' $design > t_new_design.txt
        awk '{gsub(/\\]/,"_")}1' t_new_design.txt > new_design.fa
                      
        zcat $fastq_bc | wc -l  > count_fastq.txt
                
        """
        
    }
}

/*
* count fastq and bam length remove the illegal regex characters
* and make design file
*/
if (params.labels == 0){
    process 'count_bc_nolab' {
        tag 'count'
        label 'shorttime'
        publishDir params.outdir, mode:'copy'
    
        input:
        file(fastq_bc) from fastq_bc
        file(design) from design
        //file(params.condaloc)
        //file(label) from labels
    
        output:
        file 'count_fastq.txt' into bc_ch
        file "new_label.txt" into fixed_label
        file "new_design.fa" into fixed_design
        """
        #!/bin/bash
        cv=\$(which conda)
        cv1=\$(dirname "\$cv")
        cv2=\$(dirname "\$cv1")
        cv3=\${cv2}"/bin/activate"
        echo \$cv3
        source \$cv3 mpraflow_py36
        
        #CREATE LABEL FILE
        awk -F'\t' 'BEGIN {OFS = FS} NR%2==1 {print substr(\$1,2,length(\$1)),"test"}' $design > label.txt
        
        #remove illegal regex characters
        awk '{gsub(/\\[/,"_")}1' label.txt > t_new_label.txt
        awk '{gsub(/\\]/,"_")}1' t_new_label.txt > new_label.txt
        
        awk '{gsub(/\\[/,"_")}1' $design > t_new_design.txt
        awk '{gsub(/\\]/,"_")}1' t_new_design.txt > new_design.fa
        
        zcat $fastq_bc | wc -l  > count_fastq.txt
        
        """
        
    }
}


/*
* STEP 1: Align
* Process 1A: create BWA reference
*/

process 'create_BWA_ref' {
    tag "make ref"
    label 'shorttime'
    input:
    file(design) from fixed_design
    file(label) from fixed_label

    output:
    file "${design}.fai" into reference_fai
    file "${design}.bwt" into reference_bwt
    file "${design}.sa" into reference_sa
    file "${design}.pac" into reference_pac
    file "${design}.ann" into reference_ann
    file "${design}.amb" into reference_amb
    file "${design}.dict" into reference_dict

    script:
    """
    #!/bin/bash
    cv=\$(which conda)
    cv1=\$(dirname "\$cv")
    cv2=\$(dirname "\$cv1")
    cv3=\${cv2}"/bin/activate"
    echo \$cv3
    source \$cv3 mpraflow_py36

    bwa index -a bwtsw $design
    samtools faidx $design
    picard CreateSequenceDictionary REFERENCE=$design OUTPUT=$design".dict"
    """
}

/*
*CHUNKING FASTQ
*/

Channel
    .fromPath(fastq_insert)
    .splitFastq( by: params.split, file: true )
    .set{ R1_ch }

if(params.fastq_insertPE != 0){
    Channel
        .fromPath(fastq_insertPE)
        .splitFastq( by: params.split, file: true )
        .set{ R3_ch }
}



/*
*Process 1B: merge Paired end reads
*/

if(params.fastq_insertPE != 0){
    process 'PE_merge' {
        tag 'merge'
        label 'shorttime'

        input:
        file(fastq_insert) from R1_ch
        file(fastq_insertPE) from R3_ch
        output:
        file "*merged.fastqjoin" into mergedPE

        script:
        """
        #!/bin/bash
        cv=\$(which conda)
        cv1=\$(dirname "\$cv")
        cv2=\$(dirname "\$cv1")
        cv3=\${cv2}"/bin/activate"
        echo \$cv3
        source \$cv3 mpraflow_py36        
        
        fastq-join $fastq_insert $fastq_insertPE -o ${fastq_insert}_merged.fastq
        
        """
       
    }
}
        

/*
* Process 1C: align with BWA
*/

//paired ends
if(params.fastq_insertPE != 0){
    process 'align_BWA_PE' {
        tag "align"
        label 'longtime'
    
        input:
        file(design) from fixed_design
        file(chunk) from mergedPE
        file(reference_fai) from reference_fai
        file reference_bwt from reference_bwt
        file reference_sa from reference_sa
        file reference_pac from reference_pac
        file reference_ann from reference_ann
        file reference_amb from reference_amb
        file reference_dict from reference_dict
    
        output:
        file "${params.out}.${chunk}.sorted.bam" into s_bam
        file '*count_bam.txt' into bam_ch
    
        script:
        """
        #!/bin/bash
        cv=\$(which conda)
        cv1=\$(dirname "\$cv")
        cv2=\$(dirname "\$cv1")
        cv3=\${cv2}"/bin/activate"
        echo \$cv3
        source \$cv3 mpraflow_py36
        
        bwa mem $design $chunk | samtools sort - -o ${params.out}.${chunk}.sorted.bam
        echo 'bam made'
        samtools view ${params.out}.${chunk}.sorted.bam | head
        samtools view ${params.out}.${chunk}.sorted.bam | wc -l > ${chunk}.count_bam.txt       

        """
    }
}

//single end
if(params.fastq_insertPE == 0){
    process 'align_BWA_S' {
        tag "align"
        label 'longtime'
    
        input:
        file(design) from fixed_design
        file(chunk) from R1_ch
        file(params.out)
        file(reference_fai) from reference_fai
        file reference_bwt from reference_bwt
        file reference_sa from reference_sa
        file reference_pac from reference_pac
        file reference_ann from reference_ann
        file reference_amb from reference_amb
        file reference_dict from reference_dict
    
        output:
        file "${params.out}.${chunk}.sorted.bam" into s_bam
        file '*${chunk}_count_bam.txt' into bam_ch
    
        script:
        """
        #!/bin/bash
        cv=\$(which conda)
        cv1=\$(dirname "\$cv")
        cv2=\$(dirname "\$cv1")
        cv3=\${cv2}"/bin/activate"
        echo \$cv3
        source \$cv3 mpraflow_py36
 
        bwa mem $design $chunk | samtools sort - -o ${params.out}.${chunk}.sorted.bam
        echo 'bam made'
        samtools view ${params.out}.${chunk}.sorted.bam | head
        samtools view ${params.out}.${chunk}.sorted.bam | wc -l > ${chunk}_count_bam.txt
 
        """
    }
}

/*
*COLLCT FASTQ CHUNCKS
*/

process 'collect_chunks'{
    label 'shorttime'
    
    input:
    file sbam_list from s_bam.collect()
    file count_bam from bam_ch.collect()       

    output:
    file 's_merged.bam' into s_merge
    file 'count_merged.txt' into ch_merge

    script:
    """
    #!/bin/bash
    cv=\$(which conda)
    cv1=\$(dirname "\$cv")
    cv2=\$(dirname "\$cv1")
    cv3=\${cv2}"/bin/activate"
    echo \$cv3
    source \$cv3 mpraflow_py36

    #collect sorted bams into one file
    samtools merge all.bam ${sbam_list} 
    samtools sort all.bam -o s_merged.bam    

    #collect bam counts into one file
   
    samtools view s_merged.bam | wc -l > count_merged.txt

    """
}




/*
* Assign barcodes to element sequences
*/

process 'map_element_barcodes' {
    tag "assign"
    label "shorttime"
    publishDir params.outdir, mode:'copy'

    input:
    params.mapq
    params.baseq
    params.cigar
    file(fastq_bc) from fastq_bc 
    file count_fastq from bc_ch
    file count_bam from ch_merge     
    file bam from s_merge

    output:
    file "${params.out}_coords_to_barcodes.pickle" into map_ch
    file "${params.out}_barcodes_per_candidate-no_repeats-no_jackpots.feather" into count_table_ch
    file "${params.out}_barcode_counts.pickle"

    script:
    """
    #!/bin/bash 
    cv=\$(which conda)
    cv1=\$(dirname "\$cv")
    cv2=\$(dirname "\$cv1")
    cv3=\${cv2}"/bin/activate"
    echo \$cv3
    source \$cv3 mpraflow_py36

    echo "test assign inputs"
    echo ${params.mapq}
    echo ${params.baseq}
    echo ${params.baseq}
    echo $fastq_bc
    zcat $fastq_bc | head
    
    echo ${count_fastq} 
    echo ${count_bam}    
    cat ${count_fastq}     
    cat ${count_bam} 
    
    python ${"$baseDir"}/src/nf_ori_map_barcodes.py ${"$baseDir"} ${fastq_bc} ${count_fastq} $bam ${count_bam} ${params.out} ${params.mapq} ${params.baseq} ${params.cigar}
    """
    

}

/*
* Filter barcodes for minimum coverage and unique mapping
*/

process 'filter_barcodes' {
    tag "$filter"
    label "shorttime"
    publishDir params.outdir, mode:'copy'
    input:
    params.min_cov
    params.out
    file(map) from map_ch
    file(table) from count_table_ch
    file(label) from fixed_label 
    output:
    file "${params.out}_filtered_coords_to_barcodes.pickle"
    file "${params.out}_original_counts.png"
    file "original_count_summary.txt"
    file "${params.out}_filtered_counts.png"
    file "filtered_count_summary.txt" 
 
    script:
    """
    #!/bin/bash 
    cv=\$(which conda)
    cv1=\$(dirname "\$cv")
    cv2=\$(dirname "\$cv1")
    cv3=\${cv2}"/bin/activate"
    echo \$cv3
    source \$cv3 mpraflow_py36   
 
    python ${"$baseDir"}/src/nf_filter_barcodes.py ${params.out} ${map} ${table} ${params.min_cov} ${params.min_frac} $label
    """


}





/*
 * Completion e-mail notification
 */
/*
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[MPRAflow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[MPRAflow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", attach1: "$baseDir/results/Documentation/pipeline_report.html", attach2: "$baseDir/results/pipeline_info//MPRAflow_report.html", attach3: "$baseDir/results/pipeline_info//MPRAflow_timeline.html" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[/MPRAflow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[/MPRAflow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[/MPRAflow] Pipeline Complete"

}

*/
