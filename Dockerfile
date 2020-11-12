FROM rstudio/r-base:4.0.2-xenial AS builder
MAINTAINER Marcin Kierczak <marcin.kierczak_ANTISPAM_scilifelab.se>

ARG GITHUB_PAT

ENV RENV_VERSION 0.12.0
ENV GITHUB_PAT=$GITHUB_PAT

RUN apt update -y && apt install -y \
libssl-dev \
libxml2-dev

RUN echo $GITHUB_PAT > ~/.Renviron
RUN mkdir project
#RUN mkdir project/scripts
#RUN mkdir project/scripts/src
#RUN mkdir project/scripts/doc
#RUN mkdir project/scripts/R

WORKDIR project/
  COPY renv.lock renv.lock

RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
RUN R -e 'renv::restore()'
RUN R -e "remotes::install_github('NBISweden/a_johansson_2020')"

#WORKDIR project/scripts/
#COPY src/* src/
#COPY doc/* doc/
#COPY R/* R/
