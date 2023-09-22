# gene-fusions-nf
A nextflow pipeline to call somatic rna fusions from RNA-seq data using arriba

## Dependencies
1. Nextflow : for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.
2. [arriba](https://github.com/suhrig/arriba) 

**A conda receipe, and docker and singularity containers are available with all the tools needed to run the pipeline (see "Usage")**

## Input
  | Type      | Description     |
  |-----------|---------------|
  | reads      | Path to input data |
  |bams        | Path to folder containing BAM files for processing |

## Parameters

  * #### Mandatory
| Name      | Example value | Description     |
|-----------|---------------|-----------------|
|ref_fa | ref.fa | Path to fasta reference |
|ref_gtf | annot.gtf | Path to GTF annotation|
|star_index | index/ | Path to STAR-Index reference|


  * #### Optional
| Name      | Default value | Description     |
|-----------|---------------|-----------------|
|reads_csv | NULL | file with tabular data for each sample to process |
|svs | NULL | file with tabular data for each sample with structural variants in bedpe format |

  * #### Flags

Flags are special parameters without value.

| Name      | Description     |
|-----------|-----------------|
|--arriba_plot | by default plot all the gene fusions detected by arriba. set to false to deactivate |
| --help    | Display help |

## Usage
  ```
  nextflow run  iarc/nf-gene-fusions -r v1.1 -profile singularity --reads 'test_dataset/reads/*.R{1,2}.fastq.gz' --ref_fa test_dataset/genome.fa --ref_gtf test_dataset/genome.gtf
  ```
      or
  ```
  nextflow run  iarc/nf-gene-fusions -r v1.1 -profile singularity --bams '/path/to/bams/' --ref_fa test_dataset/genome.fa --ref_gtf test_dataset/genome.gtf
  ```
  
To run the pipeline without singularity just remove "-profile singularity". Alternatively, one can run the pipeline using a docker container (-profile docker) the conda receipe containing all required dependencies (-profile conda).

## Output
  | Type      | Description     |
  |-----------|---------------|
  | outdir   | Folder with fusion genes file |



<!--- ## Detailed description (optional section) --->

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Alex Di Genova    |   | Developer |

## References
Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research. March 2021 31: 448-460; Published in Advance January 13, 2021. doi: 10.1101/gr.257246.119


