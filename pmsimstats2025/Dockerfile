# syntax=docker/dockerfile:1.4
#======================================================================
# ZZCOLLAB Ubuntu X11 Minimal Profile
#======================================================================
# Profile: ubuntu_x11_minimal (~1.0GB)
# Purpose: Minimal R with X11 support for graphics viewing
# Base: rocker/r-ver (R + minimal dependencies)
# Packages: renv (binary from r2u)
#
# Build: DOCKER_BUILDKIT=1 docker build \
#          -f Dockerfile.ubuntu_x11_minimal \
#          -t myteam/project:x11-minimal .
# Run: Requires XQuartz (macOS) or X11 (Linux)
#======================================================================

# ARGs before FROM are available only for FROM instruction
ARG R_VERSION=4.5.1

FROM rocker/r-ver:4.5.1

# Re-declare ARGs after FROM for use in build stages
ARG USERNAME=analyst
ARG RENV_VERSION=1.1.5
ARG DEBIAN_FRONTEND=noninteractive

# Reproducibility-critical environment variables (Pillar 1)
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    TZ=UTC \
    OMP_NUM_THREADS=1 \
    RENV_PATHS_CACHE=/home/analyst/.cache/R/renv

# Install runtime dependencies for X11
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt/lists,sharing=locked \
    set -ex && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
        xauth \
        libx11-6 \
        libxt6 \
        libcairo2 \
        libfontconfig1 \
        libfreetype6 \
        libpng16-16 \
        libjpeg8

# Configure R to use CRAN
RUN echo "options(repos = c(CRAN = 'https://cloud.r-project.org'))" \
        >> /usr/local/lib/R/etc/Rprofile.site

# Install renv from CRAN
RUN --mount=type=cache,target=/tmp/R-cache/4.5.1 \
    R -e "install.packages('renv')"

# Create non-root user and set up environment
RUN useradd --create-home --shell /bin/bash analyst && \
    chown -R analyst:analyst /usr/local/lib/R/site-library && \
    mkdir -p /home/analyst/.cache/R/renv && \
    chown -R analyst:analyst /home/analyst/.cache

# Switch to non-root user
USER analyst

# WORKDIR automatically creates directory with correct ownership
WORKDIR /home/analyst/project

# Health check: Verify R is functional
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=2 \
    CMD R --quiet --slave -e "quit(status = 0)" || exit 1

# Note: Project files mounted at runtime via -v $(pwd):/home/analyst/project
# This keeps image reusable across projects. Run renv::restore() in first session.
# For project-specific images, add COPY renv.lock and RUN renv::restore() before USER

CMD ["R", "--quiet"]
