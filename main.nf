#!/usr/bin/env nextflow

//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarcbioinfo/gene-fusions-nf --reads '*_R{1,2}.fastq.gz'  -profile singularity

    Mandatory arguments:
      --reads [file]                Path to input data
      --bams  [folder]              Path to folder containing BAM files for processing

    References
      --ref_fa [file]                 Path to fasta reference
      --ref_gtf [file]                Path to GTF annotation
      --star_index [file]             Path to STAR-Index reference

    Input alternatives:
      --reads_csv                   file with tabular data for each sample to process [sampleID fwd_path rev_path]
      --svs                         file with tabular data for each sample with structural variants in bedpe format [sampleID sv_path]

    Visualization :

    --arriba_plot [bool]          by default plot all the gene fusions detected by arriba. set to false to inactivate.

      Test dataset:

      The subdirectory test_dataset contains a small simulated dataset to test the whole workflow.

      Run examples:

      nextflow run  iarc/nf-gene-fusions -r v1.1 -singularity --reads 'test_dataset/reads/*.R{1,2}.fastq.gz' --ref_fa test_dataset/genome.fa --ref_gtf test_dataset/genome.gtf
      or
      nextflow run  iarc/nf-gene-fusions -r v1.1 -singularity --bams '/path/to/bams/' --ref_fa test_dataset/genome.fa --ref_gtf test_dataset/genome.gtf

      """.stripIndent()
}


if(params.bams){mode="BAM"}else{mode="reads"}

process build_star_index {
    tag "star-index"
    label 'load_medium'

    publishDir params.outdir, mode: 'copy'

    input:
        path(fasta) 
        path(gtf) 

    output:
        path("star-index"), emit: index

    when: !(params.star_index)

    script:
     //adjust a variable for working with a smaller reference
    def opt_test = params.debug ? "--genomeSAindexNbases 8" : ''
    if(params.debug){
      """
      mkdir star-index
      """
    }else{
    """
    mkdir star-index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${params.read_length - 1} \\
        --genomeDir star-index/ \\
        --genomeFastaFiles ${fasta} \\
        ${opt_test}
    """
  }
}


//we map the reads only if the bams are not mapped
process star_mapping{
  tag "${sample}"
  label 'load_low2'
  //we can remove this to avoid keeping the bam files
  publishDir "${params.outdir}/star_mapping", mode: 'copy', pattern: "${sample}.{Log.final.out,ReadsPerGene.out.tab}"

  input:
      tuple val(sample), file(reads) 
      file(star_index) 
  output:
      //star bam files
      tuple val(sample), file("${sample}_STAR.bam"), emit: bams
      //star mapping stats and gene counts *.{tsv,txt}
      tuple val(sample), file("${sample}.{Log.final.out,ReadsPerGene.out.tab,fastq.log}") optional true

  script:
  if(params.debug){
    """
    touch ${sample}_STAR.bam
    touch ${sample}.Log.final.out
    touch ${sample}.ReadsPerGene.out.tab
    """
  //BAMs are given as input
  }else if(mode=="BAM"){
  """
  #we map the reads by picking pairs from input BAM
  samtools collate -un 10 -o paired_${sample}.bam ${reads} tmp_${sample}
  samtools fastq -1 paired1.fq.gz -2 paired2.fq.gz -0 /dev/null -s /dev/null paired_${sample}.bam > ${sample}.fastq.log
  rm -f paired_${sample}.bam

  STAR \\
   --runThreadN ${task.cpus} \\
   --genomeDir ${star_index} \\
   --genomeLoad NoSharedMemory \\
   --readFilesIn  paired1.fq.gz paired2.fq.gz\\
   --readFilesCommand zcat \\
   --outSAMtype BAM Unsorted --outSAMunmapped Within \\
   --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \\
   --chimSegmentMin 10 --chimOutType WithinBAM SoftClip \\
   --chimJunctionOverhangMin 10 \\
   --chimScoreMin 1 --chimScoreDropMax 30 \\
   --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \\
   --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 \\
   --quantMode GeneCounts \\
   --outFileNamePrefix ${sample}.

  #we rename the defaul star output
  mv ${sample}.Aligned.out.bam ${sample}_STAR.bam
  #we delete the paired reads
  rm -f paired1.fq.gz paired2.fq.gz
  """
  //normal reads are given as input
  }else{
  """
  STAR \\
   --runThreadN ${task.cpus} \\
   --genomeDir ${star_index} \\
   --genomeLoad NoSharedMemory \\
   --readFilesIn ${reads} \\
   --readFilesCommand zcat \\
   --outSAMtype BAM Unsorted --outSAMunmapped Within \\
   --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \\
   --chimSegmentMin 10 --chimOutType WithinBAM SoftClip \\
   --chimJunctionOverhangMin 10 \\
   --chimScoreMin 1 --chimScoreDropMax 30 \\
   --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \\
   --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 \\
   --quantMode GeneCounts \\
   --outFileNamePrefix ${sample}.

  #we rename the defaul star output
  mv ${sample}.Aligned.out.bam ${sample}_STAR.bam
  """
 }
  /*
  #The STAR gene counts coincide with those produced by htseq-count with default parameters.
  #The file colums are:
  #column 1: gene ID
  #column 2: counts for unstranded RNA-seq
  #column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse
  */
}
//we create a channel with bam and sv files

//}
/*
 * run arriba fusion with genomic SVs
 * In case of the Variant Call Format, the file must comply with the VCF specification for structural variants.
 * In particular, Arriba requires that the SVTYPE field be present in the INFO column and specify one of the four values BND, DEL, DUP, INV.
 * In addition, for all SVTYPEs other than BND, the END field must be present and specify the second breakpoint of the structural variant.
 * Structural variants with single breakends are silently ignored.
*/
process arriba_sv {
    tag "${sample}-arriba_sv"
    label 'load_low1'

    publishDir "${params.outdir}/arriba/", mode: 'copy'

    input:
        tuple val(sample), file(vcf), file(bam)
        //path(arriba_lib)
        file(fasta)
        file(gtf)

    output:
        tuple val(sample), file("${sample}_arriba_sv.tsv") optional true
        file("*.{tsv,txt,log}")
    //we use this methods when SVs are given

    script:
    def extra_params = params.arriba_opt ? params.arriba_opt : ''
    if(params.debug){
      """
      echo arriba \\
          -x ${bam} \\
          -a ${fasta} \\
          -g ${gtf} \\
          -d ${vcf} \\
          -b ${params.arriba_lib}/blacklist_hg38_GRCh38_v2.1.0.tsv.gz \\
          -o ${sample}_arriba_sv.tsv -O ${sample}_discarded_arriba_sv.tsv \\
          ${extra_params}
          touch ${sample}_arriba_sv.tsv
      """
    }else{
    """
    arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        -d ${vcf} \\
        -b ${params.arriba_lib}/blacklist_hg38_GRCh38_v2.1.0.tsv.gz \\
        -o ${sample}_arriba_sv.tsv -O ${sample}_discarded_arriba_sv.tsv \\
        ${extra_params}  > ${sample}_arriba.log
    """
  }
}



/*
 * run arriba fusion without genomic structural variants
*/

process arriba {
    tag "${sample}-arriba"
    label 'load_low1'

    publishDir "${params.outdir}/arriba/", mode: 'copy'

    input:
        tuple val(sample), file(bam) 
        //path(arriba_lib) 
        file(fasta) 
        file(gtf) 

    output:
        tuple val(sample), file("${sample}_arriba.tsv") optional true 
        file("*.{tsv,txt,log}") 
    //we use this methods when no SVs are given

    script:
    def extra_params = params.arriba_opt ? params.arriba_opt : ''

    if(params.debug){
    """
    echo arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        -b ${params.arriba_lib}/blacklist_hg38_GRCh38_v2.1.0.tsv.gz \\
        -o ${sample}_arriba.tsv -O ${sample}_discarded_arriba.tsv \\
        ${extra_params}
      touch ${sample}_arriba.tsv
    """
   }else{
     """
     arriba \\
         -x ${bam} \\
         -a ${fasta} \\
         -g ${gtf} \\
         -b ${params.arriba_lib}/blacklist_hg38_GRCh38_v2.1.0.tsv.gz \\
         -o ${sample}_arriba.tsv -O ${sample}_discarded_arriba.tsv \\
         ${extra_params}  > ${sample}_arriba.log
    """
   }
}


/*
 * Arriba plot
 *
 */
process arriba_visualization {
    tag "${sample}-plot-fusion"
    label 'load_low1'

    publishDir "${params.outdir}/arriba/plot", mode: 'copy'

    input:
        //path(arriba_lib)
        file(gtf)
        tuple val(sample), file(bam), file(fusions) 
    output:
        file("${sample}.pdf") optional true

    when: params.arriba_plot

    script:
     //we do not plot the cytobans and protein domains for the test
    def opt_test = params.debug ? "" : "--cytobands=${params.arriba_lib}/cytobands_hg38_GRCh38_v2.1.0.tsv --proteinDomains=${params.arriba_lib}/protein_domains_hg38_GRCh38_v2.1.0.gff3"
    if(params.debug){
      //bams shold be sorted because were mapped with STAR
      """
      echo samtools sort -@ ${task.cpus} -O bam ${bam} > Aligned.sortedByCoord.out.bam
      echo samtools index Aligned.sortedByCoord.out.bam
      echo draw_fusions.R \\
          --fusions=${fusions} \\
          --alignments=Aligned.sortedByCoord.out.bam \\
          --output=${sample}.pdf \\
          --annotation=${gtf} \\
          ${opt_test}
        touch ${sample}.pdf

      """
  }else{
    """
    samtools sort -@ ${task.cpus} -O bam ${bam} > Aligned.sortedByCoord.out.bam
    samtools index Aligned.sortedByCoord.out.bam
    draw_fusions.R \\
        --fusions=${fusions} \\
        --alignments=Aligned.sortedByCoord.out.bam \\
        --output=${sample}.pdf \\
        --annotation=${gtf} \\
        ${opt_test}
    #we remove the tmp files
    rm -f Aligned.sortedByCoord.out.bam
    rm -f Aligned.sortedByCoord.out.bam.bai
    """
  }
}

//main workflow
workflow{
//display help information
  if (params.help){ show_help(); exit 0;}
  //display the header of the tool
  log.info IARC_Header()
  log.info tool_header()

//check star index or fasta and GTF
if (!params.star_index && (!params.ref_fa && !params.ref_gtf)) exit 1, "Either specify a STAR-INDEX or Fasta and GTF files!"

// create input channels
ch_fasta = Channel.value(file(params.ref_fa)).ifEmpty{exit 1, "Fasta file not found: ${params.ref_fa}"}
ch_gtf = Channel.value(file(params.ref_gtf)).ifEmpty{exit 1, "GTF annotation file not found: ${params.ref_gtf}"}

//init arriba values
params.arriba_lib = "/opt/conda/envs/gene-fusions/var/lib/arriba"
//adding arriba lib
//arriba_lib = Channel.value(file(params.arriba_lib)).ifEmpty{exit 1, "Arriba lib directory not found!"}

//expect a file with header "sampleID fwd_path rev_path"
//see file ./test_dataset/sample_fwrev.txt
if(params.reads_csv) {
      read_files_star = Channel.fromPath(file(params.reads_csv)).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row.sampleID, [file(row.fwd), file(row.rev)]]}
                      .ifEmpty{exit 1, "params.reads_csv was empty - no input files supplied" }
      //read_files_star.view()

}else if(params.bams){
    //we process the bam files
      if (file(params.bams).listFiles().findAll { it.name ==~ /.*cram/ }.size() > 0){
        	println "BAM files found, proceed with arriba";
        read_files_star=Channel.fromPath( params.bams+'/*.cram' )
           .map {  path -> [ path.name.replace(".cram",""), path] }//.view()
          mode="BAM"
        }else{
        	println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
        }

}else{
    //expect a regular expresion like '*_{1,2}.fastq.gz'
    read_files_star = Channel.fromFilePairs(params.reads, size: 2 )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
}

//sv are given as input to include in arriba
if (params.svs){
      println "SV file provided";
      //Channel for vcf files
      vcf_files = Channel.fromPath(file(params.svs)).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row.sampleID, file(row.bedpe)]}
                      .ifEmpty{exit 1, "params.reads_svs was empty - no vcf files supplied" }
                      //.view()                    
}


// Build STAR index
if(mode!="BAM"){
  build_star_index(ch_fasta,ch_gtf)
  ch_star_index = params.star_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "STAR index not found: ${params.star_index}" } : build_star_index.out.index
  //ch_star_index = ch_star_index
  //map the rna-seq reads to the genome
  star_mapping(read_files_star,ch_star_index)
  read_files_star = star_mapping.out.bams
}
//run arriba
if(params.svs){
  println "Run arriba in SV mode";
  star_bam_sv = vcf_files.join(read_files_star).view()
  arriba_sv(star_bam_sv,ch_fasta,ch_gtf)
}else{
  println "Run arriba";
  arriba(read_files_star,ch_fasta,ch_gtf)
}

//we merge into a single channel the arriba result + the star mapping
plot_arriba = params.svs ? read_files_star.join(arriba_sv.out[0]):read_files_star.join(arriba.out[0])

arriba_visualization(ch_gtf,plot_arriba)
}


//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        F\u001b[31;1mU\u001b[32;1mS\u001b[33;1mI\u001b[0mO\u001b[33;1mN\u001b[31;1m : Gene\u001b[32;1m Fusion\u001b[33;1m Caller\u001b[31;1m (${workflow.manifest.version})
        """
}
