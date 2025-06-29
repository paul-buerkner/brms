FROM rocker/verse:latest

LABEL maintainer="Sermet Pekin <sermet.pekin@gmail.com>"

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libglpk-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev libgmp-dev build-essential

# Install R packages from CRAN
RUN Rscript -e "install.packages(c( \
    'remotes', 'rcmdcheck', 'testthat', 'knitr', 'rmarkdown', 'R.rsp', \
    'lme4', 'mice', 'BH', 'MCMCglmm', 'ape', 'RcppEigen', 'emmeans', \
    'projpred', 'priorsense', 'RWiener', 'rtdists', 'extraDistr', \
    'spdep', 'mnormt', 'splines2', 'arm', 'statmod', 'diffobj', 'betareg' \
))"

# Install brms from your GitHub fork (optional)
# RUN Rscript -e "remotes::install_github('SermetPekin/brms')"

# Install cmdstanr and CmdStan
RUN Rscript -e "remotes::install_github('stan-dev/cmdstanr', dependencies = TRUE)"
RUN Rscript -e "cmdstanr::install_cmdstan(cores = parallel::detectCores())"

# Set default workdir
WORKDIR /home/rstudio
