#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarc/nf-gene-fusions --reads '*_R{1,2}.fastq.gz'  -profile docker

    Mandatory arguments:
      --reads [file]                Path to input data
      --reads_csv                   file with tabular data for each sample to process [sampleID fwd_path rev_path]
      --read_svs                    file with tabular data for each sample to process including structural variants (vcf_format) [sampleID fwd_path rev_path sv_path]

      -profile [str]                Configuration profile to use. Can use multiple (comma separated)
                                    Available: docker, singularity, test

    References
    --fasta [file]                  Path to fasta reference
    --gtf [file]                    Path to GTF annotation
    --star_index [file]             Path to STAR-Index reference

      """.stripIndent()


}

// we star coding the pipeline

// Show help message
if (params.help) exit 0, show_help()

//check star index or fasta and GTF
if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Either specify a STAR-INDEX or Fasta and GTF files!"





/*
 * Channel for reads to process
 */
if(params.reads_csv) {
      //expect a file like "sampleID fwd_path rev_path"
      reads_csv = file(params.reads_csv)
      Channel.fromPath(reads_csv).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row[0], [file(row[1]), file(row[2])]]}
                      .ifEmpty{exit 1, "params.reads_csv was empty - no input files supplied" }
                      .into{read_files_arriba}
      /*
        Channel.from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty{exit 1, "params.readPaths was empty - no input files supplied" }
            .into{read_files_arriba}
      */
} else {
  //expect a regular expresion like '*_{1,2}.fastq.gz'
    Channel.fromFilePairs(params.reads, size: 2 )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .into{read_files_arriba}
}


/*
 * Build STAR index
 */

process build_star_index {
    tag "${fasta}-${gtf}"
    label 'process_medium'

    publishDir params.outdir, mode: 'copy'

    input:
        file(fasta) from ch_fasta
        file(gtf) from ch_gtf

    output:
        file("star-index") into star_index

    when: !(params.star_index)

    script:
    def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    mkdir star-index
    STAR \\
        --runMode genomeGenerate \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${params.read_length - 1} \\
        --genomeDir star-index/ \\
        --genomeFastaFiles ${fasta} \\
        ${avail_mem}
    """
}

ch_star_index = params.star_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "STAR index not found: ${params.star_index}" } : star_index

ch_star_index = ch_star_index.dump(tag:'ch_star_index')


/*
 * run arriva fusion caller
 */
process arriba {
    tag "${sample}"
    label 'process_medium'

    publishDir "${params.outdir}/Arriba/${sample}", mode: 'copy'

    input:
        set val(sample), file(reads) from read_files_arriba
        file(reference) from reference.arriba
        file(star_index) from ch_star_index
        file(fasta) from ch_fasta
        file(gtf) from ch_gtf

    output:
        set val(sample), file("${sample}_arriba.tsv") optional true into arriba_tsv
        set val(sample), file("${sample}_arriba.bam") optional true into arriba_bam
        file("*.{tsv,txt}") into arriba_output


    script:
    def extra_params = params.arriba_opt ? params.arriba_opt : ''
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --genomeLoad NoSharedMemory \\
        --readFilesIn ${reads} \\
        --readFilesCommand zcat  \\
        --outSAMtype BAM Unsorted --outSAMunmapped Within \\
        --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \\
        --chimSegmentMin 10 --chimOutType WithinBAM SoftClip \\
        --chimJunctionOverhangMin 10 --chimScoreMin 1 \\
        --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentReadGapMax 3 |
        tee star-arriba.out.bam |

    arriba \\
        -x /dev/stdin \\
        -a ${fasta} \\
        -g ${gtf} \\
        -b ${arriba_lib}/blacklist_hg38_GRCh38_v2.0.0.tsv.gz \\
        -o ${sample}_arriba.tsv -O ${sample}_discarded_arriba.tsv \\
        -T -P ${extra_params}

    mv star-arriba.out.bam ${sample}_arriba.bam
    """
}
//arriba visualization
arriba_visualization = arriba_bam.join(arriba_tsv)

/*
 * Arriba Visualization
 */
process arriba_visualization {
    tag "${sample}"
    //label 'process_medium'

    publishDir "${params.outdir}/Arriba/${sample}", mode: 'copy'

    input:
        file(reference) from reference.arriba_vis
        set sample, file(bam), file(fusions) from arriba_visualization
        file(gtf) from ch_gtf

    output:
        file("${sample}.pdf") optional true into arriba_visualization_output

    when: params.arriba_vis

    script:
    """
    samtools sort -@ ${task.cpus} -O bam ${bam} > Aligned.sortedByCoord.out.bam
    samtools index Aligned.sortedByCoord.out.bam
    draw_fusions.R \\
        --fusions=${fusions} \\
        --alignments=Aligned.sortedByCoord.out.bam \\
        --output=${sample}.pdf \\
        --annotation=${gtf} \\
        --cytobands=${arriba_lib}/cytobands_hg38_GRCh38_v2.0.0.tsv \\
        --proteinDomains=${arriba_lib}/protein_domains_hg38_GRCh38_v2.0.0.gff3
    """
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
# Nextflow pilelines for cancer genomics.########################################
"""
}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        F\u001b[31;1m \\
        U\u001b[32;1m \\
        S\u001b[33;1m \\
        I\u001b[0m \\
        O\u001b[33;1m \\
        N\u001b[31;1m : \\
        Gene\u001b[32;1m \\
        Fusion\u001b[33;1m \\
        Caller\u001b[31;1m \\
"""
}
