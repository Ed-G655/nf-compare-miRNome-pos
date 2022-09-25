
################## BASE IMAGE ######################
FROM rocker/tidyverse
################## METADATA ######################
LABEL base_image="rocker/tidyverse"
LABEL version="0"
LABEL software="nf-compare-miRNome"
LABEL software.version="0.2"
LABEL about.summary="Nextflow pipeline to extract and create miRNA and 3'UTR consensus sequences for analysis \n with targetscan and miRmap"
LABEL about.home="https://github.com/Ed-G655/nf-miRNA-SNPs-classify.nf.git"
LABEL about.documentation="https://github.com/Ed-G655/nf-miRNA-SNPs-classify.nf.git"
LABEL about.tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Jose Eduardo Garcia Lopez <jeduardogl655@gmail.com>

################## INSTALLATION ######################

# install pipeline requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
  wget \
  curl \
  gcc-multilib \
  apt-utils \
  curl \
  perl \
  zip \
  unzip

# install R packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  && install2.r --error \
    vroom \
    cowplot \
    eulerr \
    ggvenn \
    pacman

RUN apt-get install -y --no-install-recommends less \
  grep

COPY ./bin /usr/local/bin
