#!/bin/sh
set -ex
wget https://github.com/Reference-LAPACK/lapack/archive/3.9.0.tar.gz
tar -xzvf lapack-3.9.0.tar.gz
cd lapack-3.9.0 && ./configure --prefix=/usr && make && sudo make install