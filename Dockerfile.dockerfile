# Use the official R Shiny base image
FROM rocker/shiny:latest

# Install system dependencies (optional – customize as needed)
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev

# Install R packages required by your app
# Modify this line to match the packages your app uses
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'plotly', 'leaflet', 'reshape2', 'leaflet', 'shinyWidgets', 'Polychrome', 'dplyr', 'tseries', 'DT', 'collapse', 'magrittr', 'tseries', 'lmtest', 'sandwich', 'jtools', 'xts', 'ggfortify', 'lubridate', 'MTE', 'quantreg', 'meboot', 'foreach', 'doParallel', 'future', 'furrr', 'purrr', 'rqPen'), repos='https://cloud.r-project.org')"

# Copy your app into the Docker container
COPY . /srv/shiny-server/

# Expose Shiny default port
EXPOSE 8264

# Run Shiny Server
CMD ["/usr/bin/shiny-server"]
