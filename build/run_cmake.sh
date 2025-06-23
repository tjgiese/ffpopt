#!/bin/bash
set -e
set -u

PREFIX=${PWD}/../local

bash ./deletebuild.sh
cmake .. \
	-DCMAKE_INSTALL_PREFIX=${PREFIX} \
        -DBUILD_PYTHON=TRUE \
        -DPython3_EXECUTABLE=`which python3`
make install
