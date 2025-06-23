#!/bin/bash
set -e
set -u

if [ ! -e deletebuild.sh ]; then
    echo "nothing to do"
fi

for f in CMakeCache.txt  CMakeFiles  cmake_install.cmake  _deps  generated  install_manifest.txt  lapacke_mangling.h  Makefile  openblas_config.h  src *~; do
    if [ -e "${f}" ]; then
	echo "rm -fr ${f}"
	rm -fr "${f}"
    fi
done

