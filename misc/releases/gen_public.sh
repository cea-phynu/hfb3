#!/bin/bash

NOCOPY=(misc/amedee_python misc/berger_python misc/docs misc/releases/gen_private.sh misc/releases/gen_pip_private.sh go.sh src/drop.cpp src/drop.h src/link.cpp src/link.h src/field_spin_orbit_fr.cpp src/field_spin_orbit_fr.h src/field_tensor_fr.cpp src/field_tensor_fr.h gen_release.sh gen_pip_release.sh misc/studies/180Hg_asym_drop misc/pip)

echo "===== Cleaning project ====="
make clean

VERSION=$(git describe --always --dirty --tags)
DIR=release-$VERSION

echo "===== Creating release in $DIR ====="

shopt -s extglob

rm -rf $DIR
mkdir $DIR
cp -r !($DIR) $DIR/

cp .zenodo.json $DIR/

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
      mv $name $i
    else
      rm $name
    fi
  fi
done

cd ..
