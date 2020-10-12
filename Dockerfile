################# BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="nf-gene-fusions"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for nf-gene-fusions"
LABEL about.home="https://github.com/adigenova/nf-gene-fusions"
LABEL about.documentation="https://github.com/adigenova/nf-gene-fusions/README.md"
LABEL about.license_file="https://github.com/adigenova/nf-gene-fusions/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER Alex Di Genova <digenova@gmail.com>

################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n nf-gene-fusions -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-gene-fusions/bin:$PATH

# install dependencies
#RUN export DEBIAN_FRONTEND=noninteractive && \
#apt-get update -y && \
#apt-get install -y --no-install-recommends build-essential samtools STAR r-base rna-star wget ca-certificates libcurl4-openssl-dev libxml2-dev && \
#Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org"); source("https://bioconductor.org/biocLite.R"); biocLite(c("GenomicRanges", "GenomicAlignments"))'

# download and compile arriba v1.2
#RUN wget -q -O - "https://github.com/suhrig/arriba/releases/download/v1.2.0/arriba_v1.2.0.tar.gz" | tar -xzf - && \
#cd arriba_v1.2.0 && make all && cp arriba /usr/bin/
#STAR version=2.7.4a
#samtools 1.9
#REF=/data/scratch/digenovaa/mesomics/rna-seq/STAR-arriba/reference/GRCh38.primary_assembly.genome.fa
#GTF=/data/scratch/digenovaa/mesomics/rna-seq/alignment-free/gencode.v34.annotation.gtf.gz
#ARRIBA=/data/scratch/digenovaa/mesomics/rna-seq/STAR-arriba/arriba_v1.2.0/arriba
#BLACKLIST=/data/scratch/digenovaa/mesomics/rna-seq/STAR-arriba/arriba_v1.2.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz
