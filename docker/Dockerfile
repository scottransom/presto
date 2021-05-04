FROM ubuntu:18.04

# Install prerequisites
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get -y install software-properties-common && \
    add-apt-repository -y universe && \
    add-apt-repository -y multiverse && \
    apt-get -y install \
    build-essential \
    software-properties-common \
    wget \
    python3 \
    python3-dev \
    python \
    libpython3.8-dev \
    python3-pip \
    python3-numpy \
    python3-six \
    python3-scipy \
    python3-astropy \
    python3-setuptools \
    gfortran \
    libglib2.0-dev \
    libpng-dev \
    curl \
    autoconf \
    autotools-dev \
    automake \
    autogen \
    libtool \
    pkg-config \
    cmake \
    csh \
    gcc-5 \
    git \
    cvs \
    libcfitsio-dev \
    pgplot5 \
    hwloc \
    python-dev \
    libfftw3-3 \
    libfftw3-bin \
    libfftw3-dev \
    libfftw3-single3 \
    libx11-dev \
    libpnglite-dev \
    libhdf5-dev \
    libhdf5-serial-dev \
    libxml2 \
    libxml2-dev \
    libltdl-dev \
    gsl-bin \
    libgsl-dev 

RUN curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py && \
    pip install future

RUN pip3 install numpy \
    scipy \
    astropy

# Install presto python scripts
ENV PRESTO /code/presto
ENV LD_LIBRARY_PATH /code/presto/lib
ADD . /code/presto

WORKDIR /code/presto/src
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
    
# Install tempo
RUN git clone https://github.com/nanograv/tempo.git && \
    cd tempo && \
    ./prepare && \
    ./configure && \
    make && \
    make install
ENV TEMPO /home/soft/tempo

# Install tempo2
WORKDIR /home/soft/tempo2
ENV TEMPO2 /home/soft/tempo2
RUN wget https://sourceforge.net/projects/tempo2/files/latest/download && \
    tar -xzvf download && \
    rm download && \
    cp -r tempo2*/T2runtime/ $TEMPO2

# Install pgplot
WORKDIR /usr/local/src
RUN wget ftp://ftp.astro.caltech.edu/pub/pgplot/pgplot5.2.tar.gz && \
    tar zxvf pgplot5.2.tar.gz && \
    rm pgplot5.2.tar.gz
WORKDIR /usr/local/pgplot
ADD docker/drivers.list /usr/local/pgplot
RUN /usr/local/src/pgplot/makemake /usr/local/src/pgplot linux g77_gcc_aout && \
    sed -i 's/FCOMPL=g77/FCOMPL=gfortran/' makefile && \
    make && \
    make cpg && \
    make clean
ENV PGPLOT_DIR=/usr/local/pgplot
ENV PGPLOT_DEV=/Xserve
    
# Install presto
WORKDIR /code/presto/src
RUN make makewisdom && \
    make prep && \
    make -j 1 && \
    make clean
ENV PATH="/code/presto/bin/:${PATH}"

