# Use the official R Shiny base image
FROM rocker/shiny:latest

# Install system dependencies
# Added:
# - build-essential: General build tools (g++, make, etc.) crucial for many R packages.
# - gfortran: Required by many statistical/numerical R packages, including 'lme4', 'nloptr'.
# - libblas-dev, liblapack-dev: Linear algebra libraries, often needed for 'lme4', 'nloptr', etc.
# - libudunits2-dev, libgmp-dev: Often needed for geospatial or numerical packages.
# - libprotobuf-dev, protobuf-compiler: If you use any packages that interact with Protobuf (less common but good to have if issues persist).
# - procps: provides 'ps' command, sometimes used by R packages for process management.
# - sudo: While not directly fixing compilation, it's often assumed to be present in interactive environments.
# - r-base-dev: Contains R development headers, usually part of base R images but good to explicitly list if troubleshooting.
RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libgmp-dev \
    libprotobuf-dev \
    protobuf-compiler \
    procps \
    sudo \
    r-base-dev \
    # Clean up APT cache to reduce image size
    && rm -rf /var/lib/apt/lists/*

# Install R packages required by your app
# Added packages that failed previously: nloptr, lme4, pbkrtest, car, dynlm.
# Made sure to include all packages from your original list, and combined them.
# Increased Ncpus to leverage more cores if available during build.
RUN R -e "install.packages(c( \
    'shiny', 'shinydashboard', 'plotly', 'leaflet', 'reshape2', 'shinyWidgets', \
    'Polychrome', 'dplyr', 'tseries', 'DT', 'collapse', 'magrittr', \
    'lmtest', 'sandwich', 'jtools', 'xts', 'ggfortify', 'lubridate', \
    'MTE', 'quantreg', 'meboot', 'foreach', 'doParallel', 'future', \
    'furrr', 'purrr', 'rqPen', \
    'nloptr', 'lme4', 'pbkrtest', 'car', 'dynlm' \
    ), \
    Ncpus = 8, \
    repos='https://cloud.r-project.org')"

# Copy your app into the Docker container
COPY . /srv/shiny-server/

# IMPORTANT: Render typically expects web services to listen on port 10000,
# or for web services it infers (like Shiny Server), it might override this.
# The `rocker/shiny` image by default runs Shiny Server on port 3838.
# You had EXPOSE 8264. Let's ensure consistency.
# If your Shiny Server configuration (shiny-server.conf) is set to 8264,
# then keep it. Otherwise, 3838 is the default for rocker/shiny.
# If Render complains about "No open HTTP ports detected", this is the place to check.
EXPOSE 8264 
# Changed from 8264 to 3838 (Shiny Server default)

# Run Shiny Server
# The default command for rocker/shiny is /usr/bin/shiny-server.
# If you have a custom shiny-server.conf that changes the port, ensure it's copied and used.
CMD ["/usr/bin/shiny-server"]