################# BASE IMAGE #####################
#FROM continuumio/miniconda3:4.10.3
FROM mambaorg/micromamba:0.15.3
##site to test docker configuration files
# https://labs.play-with-docker.com/
################## METADATA #######################
LABEL base_image="mambaorg/micromamba"
LABEL version="0.15.3"
LABEL software="nf-gene-fusions"
LABEL software.version="1.1"
LABEL about.summary="Container image containing all requirements for nf-gene-fusions"
LABEL about.home="https://github.com/adigenova/nf-gene-fusions"
LABEL about.documentation="https://github.com/adigenova/nf-gene-fusions/README.md"
LABEL about.license_file="https://github.com/adigenova/nf-gene-fusions/LICENSE.txt"
LABEL about.license="GNU-3.0"

################## MAINTAINER ######################
MAINTAINER Alex Di Genova <digenova@gmail.com>
################## INSTALLATION ######################
USER root
#the next command install the ps command needed by nexflow to collect run metrics
RUN apt-get update && apt-get install -y procps
USER micromamba
COPY --chown=micromamba:micromamba environment.yml /tmp/environment.yml
RUN micromamba create -y -n gene-fusions -f /tmp/environment.yml && \
    micromamba clean --all --yes
ENV PATH /opt/conda/envs/gene-fusions/bin:$PATH
#arriba lib
#backlisted hg38 regions
#/opt/conda/envs/gene-fusions/var/lib/arriba/blacklist_hg38_GRCh38_v2.1.0.tsv.gz
#protein domains
#/opt/conda/envs/gene-fusions/var/lib/arriba/protein_domains_hg38_GRCh38_v2.1.0.gff3
