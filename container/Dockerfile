FROM mambaorg/micromamba:1.5.1

RUN ldd --version | head -n 1

LABEL authors="KFC" \
      description="Docker image containing R packages for spatial-ATAC-RNA analysis"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    procps \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install dependencies using micromamba
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Add R-specific system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up R environment variables
ENV PATH="/opt/conda/bin:$PATH"

# Install additional R packages that might not be available in conda
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='http://cran.rstudio.com/'); \
    BiocManager::install(c('GenomeInfoDb')); \
    if (!requireNamespace('Signac', quietly = TRUE)) remotes::install_github('stuart-lab/signac'); \
    if (!requireNamespace('ArchR', quietly = TRUE)) remotes::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories());"

# Set working directory
WORKDIR /data

# Set default command
CMD ["/bin/bash"]
