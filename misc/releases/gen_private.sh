#!/bin/bash

NOCOPY=(misc/releases)

echo "=====  Cleaning project ====="
make clean

VERSION=$(git describe --always --dirty --tags)
DIR=release-$VERSION

echo "===== Creating release in $DIR ====="

shopt -s extglob

rm -rf $DIR
mkdir $DIR
cp -r !($DIR) $DIR/

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
    cat $i | grep -I -v '#NEVER_IN_RELEASE' | sed 's/1.0.0rc21/'$VERSION'/g' > $name
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


