# ─────────────────────────────────────────────────────────────────────────────
# bioBakery Visualization Workflow — Docker image
#
# Extends the official biobakery/workflows image with any extra packages
# needed for PDF generation and report rendering.
#
# Build:
#   docker build -t biobakery/workflows-vis:latest environments/
#
# Push to registry (e.g., Docker Hub or ECR):
#   docker tag biobakery/workflows-vis:latest <registry>/<repo>:<tag>
#   docker push <registry>/<repo>:<tag>
#
# Then set in nextflow.config:
#   params.container = '<registry>/<repo>:<tag>'
# ─────────────────────────────────────────────────────────────────────────────

# Official biobakery image — includes Python, R, anadama2, biobakery_workflows,
# humann, metaphlan, kneaddata, and their R/Python dependencies.
FROM biobakery/workflows:latest

LABEL maintainer="bioBakery workflows development team <biobakery-users@googlegroups.com>"
LABEL description="bioBakery visualization workflow container (Nextflow edition)"
LABEL version="4.0.0"

USER root

# ── System dependencies ───────────────────────────────────────────────────────
# pandoc  : Pweave → PDF/HTML conversion
# texlive : LaTeX PDF backend (used by Pandoc)
# fonts   : required by matplotlib / LaTeX for clean plot rendering
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        pandoc \
        texlive-xetex \
        texlive-latex-extra \
        texlive-fonts-recommended \
        texlive-fonts-extra \
        lmodern \
        fonts-dejavu \
        curl \
        procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── Python packages ───────────────────────────────────────────────────────────
# pweave       : template engine used by biobakery_workflows document generation
# hclust2      : hierarchical clustering heatmaps (called inside Pweave templates)
# scipy ≥ 1.7  : PCoA ordination (required by biobakery_workflows/visualizations.py)
RUN pip install --no-cache-dir \
        "pweave>=0.30.3" \
        "hclust2>=1.0.0" \
        "scipy>=1.7.0" \
        "matplotlib>=3.3.0" \
        "pandas>=1.2.0" \
        "numpy>=1.20.0"

# ── R packages ────────────────────────────────────────────────────────────────
# These are required by alpha_diversity.R and beta_diversity.R.
# The base biobakery image may already include them; this ensures they are present.
RUN Rscript -e "
    pkgs <- c('vegan','ggplot2','optparse','gridExtra','dplyr','tibble','purrr','readr')
    missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
    if (length(missing) > 0) {
        install.packages(missing, repos='https://cloud.r-project.org', quiet=TRUE)
    }
"

# ── Verify installation ───────────────────────────────────────────────────────
RUN biobakery_workflows --help 2>&1 | grep -q "workflows" || \
    (echo "ERROR: biobakery_workflows not found in image" && exit 1)

RUN pandoc --version | head -1 || (echo "ERROR: pandoc not found" && exit 1)

# ── Working directory ─────────────────────────────────────────────────────────
WORKDIR /work

# Default command (overridden by Nextflow)
CMD ["biobakery_workflows", "--help"]
