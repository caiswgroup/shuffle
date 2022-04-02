#!/bin/sh
cd kissat
./configure
make all || exit 1
# build/tissat || exit 1
exec install -s build/kissat ../../bin/
