#!/bin/sh
set -ex
wget https://github.com/jchristopherson/ferror/archive/1.3.0.tar.gz
tar -xzvf ferror-1.3.0.tar.gz
cd ferror-1.3.0 && ./configure --prefix=/usr && make && sudo make install