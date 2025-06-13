#!/bin/bash

if [ ! -z "$1" ]
then
nb=$1
else
nb=$(nproc --all)
fi

echo "installing tox (dependency)"
pip install tox
if [ $? -ne 0 ]; then
  echo "could not install tox with 'pip install tox'. Maybe use a venv ?"
  exit
fi

echo "launching Python tests using $nb jobs"
make -j$nb lib/libhfb3.a && (cd misc/python; tox)

