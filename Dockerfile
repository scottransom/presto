FROM ubuntu:20.04

# Install prerequisites
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update -qq && \
    apt-get -y --no-install-recommends install \
    autoconf \
    automake \
    build-essential \
    gfortran \
    git \
    latex2html \
    libcfitsio-bin \
    libcfitsio-dev \
    libfftw3-bin \
    libfftw3-dev \
    libglib2.0-dev \
    libpng-dev \
    libtool \
    libx11-dev \
    pgplot5 \
    python3-dev \
    python3-numpy \
    python3-pip \
    tcsh \
    wget && \
    apt-get clean all && \
    rm -r /var/lib/apt/lists/*

# Add pgplot environment variables
ENV PGPLOT_DIR=/usr/local/pgplot
ENV PGPLOT_DEV=/Xserve

# Install python dependancies
RUN pip3 install numpy \
    scipy \
    astropy

# Install presto python scripts
ENV PRESTO /code/presto
ENV LD_LIBRARY_PATH /code/presto/lib
ADD . /code/presto

WORKDIR /code/presto/src
# The following is necessary if your system isn't Ubuntu 20.04
RUN make cleaner
# Now build from scratch
RUN make libpresto slalib
WORKDIR /code/presto
RUN pip3 install /code/presto && \
    sed -i 's/env python/env python3/' /code/presto/bin/*py && \
    python3 tests/test_presto_python.py 


# Installs all the C dependancies -----------------------------
WORKDIR /home/soft

# Install psrcat
RUN wget https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz && \
    gunzip psrcat_pkg.tar.gz && \
    tar -xvf psrcat_pkg.tar && \
    rm psrcat_pkg.tar && \
    cd psrcat_tar && \
    ls && \
    bash makeit && \
    cp psrcat /usr/bin
ENV PSRCAT_FILE /home/soft/psrcat_tar/psrcat.db
    
# Install tempo
RUN git clone https://github.com/nanograv/tempo.git && \
    cd tempo && \
    ./prepare && \
    ./configure && \
    make && \
    make install
ENV TEMPO /home/soft/tempo
 
# Install presto
WORKDIR /code/presto/src
RUN make makewisdom && \
    make prep && \
    make -j 1 && \
    make clean
ENV PATH="/code/presto/bin/:${PATH}"

