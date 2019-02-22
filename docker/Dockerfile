FROM ubuntu:16.04

RUN apt-get update

RUN apt-get -qq -y install pkg-config flex bison git python gcc g++

RUN apt-get -y upgrade

RUN git clone https://github.com/alexandre-lefort/these.git

RUN bin/bash -c "cd /these/ibex-lib && \
                 ./waf configure --with-optim --with-optim-MinMax --with-affine && \
                 ./waf install && \
                 cd ../../ "

RUN bin/bash -c "cd /these/tests/sm_2_barres/ibex/ && make"

RUN bin/bash -c "cd /these/tests/sm_2_barres/ && ./sm_2_barres"

