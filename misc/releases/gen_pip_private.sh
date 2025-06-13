#!/bin/bash

NOCOPY=(misc/releases)

echo "===== Cleaning project ====="
make clean

if [ -z "$1" ]
then
VERSION=$(git describe --always --dirty --tags)
else
VERSION=$1
fi

DIR=release-pip-$VERSION

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
    cat $i | grep -v '#NEVER_IN_RELEASE' | sed 's/1.0.0rc21/'$VERSION'/g' > $name
    if ! cmp -s $name $i
    then
      echo "patching $i"
      mv -f $name $i
    else
      rm $name
    fi
  fi
done

cd ..
