#!/bin/bash

if [ -z "$1" ]
then
  echo "Usage: $0 PYTHON_VERSION [MODULE_VERSION]"
  exit 0
fi


NOCOPY=(misc/amedee_python misc/berger_python misc/docs misc/releases src/drop.cpp src/drop.h src/link.cpp src/link.h src/field_spin_orbit_fr.cpp src/field_spin_orbit_fr.h src/field_tensor_fr.cpp src/field_tensor_fr.h)

echo "===== Cleaning project ====="
make clean

PYTHON_VERSION=$1

if [ -z "$2" ]
then
VERSION=$(git describe --always --dirty --tags)
else
VERSION=$2
fi

DIR=release-pip-$VERSION-$PYTHON_VERSION

echo "===== Creating PIP release in $DIR ====="

shopt -s extglob

rm -rf $DIR
mkdir $DIR
cp -r misc/pip/* $DIR/
cp -r misc/python/hfb3.i $DIR/
cp -r src $DIR/
mkdir $DIR/hfb3
mkdir $DIR/lib
mkdir $DIR/misc
cp -r examples $DIR/hfb3/
cp -r misc/deps/ $DIR/misc/
cp Makefile    $DIR/
cp Makefile.in $DIR/

cd $DIR

echo "===== Removing files ====="
for i in ${NOCOPY[@]}
do
  echo removing $i
  rm -rf $i
done

echo "===== Patching files ====="
for i in $(find . -type f)
do
  if [ -n "$(file $i | grep text)" ]
  then
    name=$(mktemp)
    if ! cmp -s $name $i
    then
      echo "patching $i"
      mv -f $name $i
    else
      rm $name
    fi
  fi
done

docker run -it --volume $(pwd):/root/ quay.io/pypa/manylinux_2_28_x86_64 /root/go$PYTHON_VERSION.sh

cd ..

