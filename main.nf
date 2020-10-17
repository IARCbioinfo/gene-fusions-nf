#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarc/nf-gene-fusions --reads '*_R{1,2}.fastq.gz'  -profile singularity

    Mandatory arguments:
      --reads [file]                Path to input data

    References
      --fasta [file]                  Path to fasta reference
      --gtf [file]                    Path to GTF annotation
      --star_index [file]             Path to STAR-Index reference

    Input alternatives:
      --reads_csv                   file with tabular data for each sample to process [sampleID fwd_path rev_path]
      --reads_svs                   file with tabular data for each sample to process including structural variants
                                    (vcf_format) [sampleID fwd_fullpath rev_fullpath sv_path].

      -profile [str]              Configuration profile to use.
                                  Available: docker and singularity

      Test dataset:

      The subdirectory test_dataset contains a small simulated dataset to test the whole workflow.

      Test run:

      nextflow run  iarc/nf-gene-fusions --reads 'test_dataset/reads/*.R{1,2}.fastq.gz' --fasta test_dataset/genome.fa --gtf test_dataset/genome.gtf

      """.stripIndent()
}

// we star coding the pipeline

// Show help message
if (params.help) exit 0, show_help()




//check star index or fasta and GTF
if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Either specify a STAR-INDEX or Fasta and GTF files!"

ch_fasta = Channel.value(file(params.fasta)).ifEmpty{exit 1, "Fasta file not found: ${params.fasta}"}
ch_gtf = Channel.value(file(params.gtf)).ifEmpty{exit 1, "GTF annotation file not found: ${params.gtf}"}

//init arriba values
arriba = [
    lib: false,
]
//adding arriba lib
arriba.lib = Channel.value(file(params.arriba_lib)).ifEmpty{exit 1, "Arriba lib directory not found!"}

log.info IARC_Header()
log.info tool_header()

/*
 * Channel for reads to process
 */
 /*
       Channel.from(params.readPaths)
           .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
           .ifEmpty{exit 1, "params.readPaths was empty - no input files supplied" }
           .into{read_files_arriba}
 */
//expect a file like "sampleID fwd_path rev_path"
if(params.reads_csv) {
      reads_csv = file(params.reads_csv)
      Channel.fromPath(reads_csv).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row[0], [file(row[1]), file(row[2])]]}
                      .ifEmpty{exit 1, "params.reads_csv was empty - no input files supplied" }
                      .set{read_files_star}
//expect a file like "sampleID fwd_path rev_path vcf_file"

if(params.reads_csv  =~ /test_dataset/){
    params.test=true;
  }
}
/*
else if (params.reads_svs){
      reads_svs = file(params.reads_svs)
      Channel.fromPath(reads_svs).splitCsv(header: true, sep: '\t', strip: true)
                      .map{row -> [ row[0], [file(row[1]), file(row[2])], file(row[3])]}
                      .ifEmpty{exit 1, "params.reads_svs was empty - no input files supplied" }
                      .set{read_files_arriba_sv}
}*/
//expect a regular expresion like '*_{1,2}.fastq.gz'
 else  {
    Channel.fromFilePairs(params.reads, size: 2 )
        .ifEmpty{exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
        .set{read_files_star}

      if(params.reads =~ /test_dataset/){
            params.test=true;
        }
}





/*
 * Build STAR index
 */

process build_star_index {
    tag "star-index"
    label 'load_medium'

    publishDir params.outdir, mode: 'copy'

    input:
        file(fasta) from ch_fasta
        file(gtf) from ch_gtf

    output:
        file("star-index") into star_index

    when: !(params.star_index)

    script:
    //def avail_mem = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : '' //${avail_mem}
    def opt_test = params.test ? "--genomeSAindexNbases 8" : ''; //adjust a variable for working with a smaller reference
    //if we are running the test we shold reduce the index size
    //def genomeSAindexNbases = task.test ? 8:14
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

ch_star_index = params.star_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "STAR index not found: ${params.star_index}" } : star_index

ch_star_index = ch_star_index.dump(tag:'ch_star_index')

//map the rna-seq reads to the genome

process star_mapping{
  tag "${sample}"
  label 'load_low2'
  //we can remove this to don't keep the bam files
  publishDir "${params.outdir}/star_mapping", mode: 'copy'

  input:
      set val(sample), file(reads) from read_files_star
      file(star_index) from ch_star_index
  output:
      set val(sample), file("${sample}_STAR.bam") into star_bam //bam
      set val(sample), file("${sample}_STAR.{result.log,gene_counts.log}") optional true into star_output //star log *.{tsv,txt}

  script:
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
      --chimJunctionOverhangMin 10  \\
    	--chimScoreMin 1 --chimScoreDropMax 30 \\
      --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 \\
    	--alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 \\
      --quantMode GeneCounts \\
      --outFileNamePrefix ${sample}.

      #we rename the defaul star output
      mv ${sample}.Aligned.out.bam ${sample}_STAR.bam
      #we rescue the mapping stats
      mv ${sample}.Log.final.out ${sample}_STAR.result.log
      #we rescue the geneCounts
      mv ${sample}.ReadsPerGene.out.tab ${sample}_STAR.gene_counts.log
  """
  /*
  #The STAR gene counts coincide with those produced by htseq-count with default parameters.
  #The file colums are:
  #column 1: gene ID
  #column 2: counts for unstranded RNA-seq
  #column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  #column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse
  */

}


/*
 * run arriba fusion
*/

process arriba {
    tag "${sample}"
    label 'process_medium'

    publishDir "${params.outdir}/Arriba/", mode: 'copy'

    input:
        set sample, file(bam) from star_bam
        file(arriba_lib) from arriba.lib
        file(fasta) from ch_fasta
        file(gtf) from ch_gtf

    output:
        set val(sample), file("${sample}_arriba.tsv") optional true into arriba_tsv
        file("*.{tsv,txt,log}") into arriba_output

    script:
    def extra_params = params.arriba_opt ? params.arriba_opt : ''
    def opt_test = params.test ? "-f blacklist" : ''; //adjust a variable for working with the smaller reference
    """
    arriba \\
        -x ${bam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        -b ${arriba_lib}/blacklist_hg38_GRCh38_v2.0.0.tsv.gz \\
        -o ${sample}_arriba.tsv -O ${sample}_discarded_arriba.tsv \\
        ${extra_params} ${opt_test} > ${sample}_arriba.log
    """
}
//arriba visualization
//arriba_visualization = arriba_bam.join(arriba_tsv)


/*
 * run arriba fusion with genomic SVs
 * In case of the Variant Call Format, the file must comply with the VCF specification for structural variants.
 * In particular, Arriba requires that the SVTYPE field be present in the INFO column and specify one of the four values BND, DEL, DUP, INV.
 * In addition, for all SVTYPEs other than BND, the END field must be present and specify the second breakpoint of the structural variant.
 * Structural variants with single breakends are silently ignored.


process arriba_sv {
    tag "${sample}"
    //label 'process_medium'

    publishDir "${params.outdir}/Arriba/${sample}", mode: 'copy'

    input:
        set val(sample), file(reads), file(vcf) from read_files_arriba_sv
        file(arriba_lib) from arriba.lib
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
        -d ${vcf} \\
        -b ${arriba_lib}/blacklist_hg38_GRCh38_v2.0.0.tsv.gz \\
        -o ${sample}_arriba.tsv -O ${sample}_discarded_arriba.tsv \\
        -T -P ${extra_params}

    mv star-arriba.out.bam ${sample}_arriba.bam
    """
}
//arriba visualization
arriba_visualization = arriba_bam.join(arriba_tsv)
*/


/*
 * Arriba Visualization

process arriba_visualization {
    tag "${sample}"
    //label 'process_medium'

    publishDir "${params.outdir}/Arriba/${sample}", mode: 'copy'

    input:
        file(arriba_lib) from arriba.lib
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
*/

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
        F\u001b[31;1mU\u001b[32;1mS\u001b[33;1mI\u001b[0mO\u001b[33;1mN\u001b[31;1m : Gene\u001b[32;1m Fusion\u001b[33;1m Caller\u001b[31;1m (${workflow.manifest.version})
        """
}
