#!/bin/bash

if [ ! -z "$1" ]
then
nb=$1
else
nb=$(nproc --all)
fi

echo "launching tests using $nb jobs"

make -j$nb tests && OMP_NUM_THREADS=1 misc/tests/gtest-parallel -w$nb bin/tests


