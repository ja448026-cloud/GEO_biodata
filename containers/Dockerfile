FROM rocker/r-ver:4.5.0

RUN R -q -e "install.packages(c('BiocManager','httr2','jsonlite','digest','yaml'), repos='https://cloud.r-project.org')" \
 && R -q -e "BiocManager::install(c('GEOquery','Biobase'), ask=FALSE, update=FALSE)"

WORKDIR /workflow
COPY . /workflow

CMD ["Rscript", "core/R/bootstrap_environment.R", "--profile", "core", "--check"]
