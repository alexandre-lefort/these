FROM ubuntu

RUN apt-get update

RUN apt-get -qq -y install pkg-config \
    flex \
    bison \
    git \
    python \
    gcc \
    g++

RUN apt-get -y upgrade

RUN git clone https://github.com/alexandre-lefort/these.git

RUN cd these/ibex-lib \
    ./waf configure --with-optim --with-optim-MinMax --with-affine \
    ./waf install \
    cd ../../ 

RUN cd /these/tests/sm_2_barres/ibex/ \
    make \
    ./sm_2_barres

