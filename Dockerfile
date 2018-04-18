FROM rocker/shiny

WORKDIR /usr/src/app

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite3-dev \
  libmariadbd-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssh2-1-dev \
  unixodbc-dev \
  && R -e "source('https://bioconductor.org/biocLite.R')" \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    dplyr \
    ggplot2 \
    devtools \
    formatR \
    remotes \
    selectr \
    caTools; apt-get update; apt-get install -y vim; apt-get install -y libncurses5-dev; apt-get install -y libbz2-dev; apt-get install -y liblzma-dev; apt-get update; apt-get install -y build-essential; apt-get install -y git; apt-get install -y gcc; apt-get install -y zlib1g-dev; apt-get install -y make;  apt-get install -y wget; git clone https://github.com/lh3/bwa.git; cd bwa;  make; wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2; tar -vxjf samtools-1.8.tar.bz2; cd samtools-1.8; ./configure; make; make install; apt-get install -y r-base; R -e "source('https://bioconductor.org/biocLite.R')"; R -e "install.packages('plotrix')"; R -e "install.packages('tidyverse')";

##how to start this worthy of being called a very bad word image
## docker run --rm -it --mount type=bind,source="$(pwd)"/8678/,target=/usr/src/app/data/ f448239880fa

CMD ["/bin/bash"]
