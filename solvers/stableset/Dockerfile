FROM ubuntu:22.04

RUN apt-get update
RUN apt install build-essential -y
RUN apt-get install openmpi-bin openmpi-doc libopenmpi-dev -y 


WORKDIR /solver

COPY . .

RUN cd external_libs/ && ./install.sh && cd /solver

RUN chmod -R 777 .

USER 1001:1001

ENV LD_LIBRARY_PATH=/tmp/mosek/8/tools/platform/linux64x86/bin

RUN make clean
RUN make

