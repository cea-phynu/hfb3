#!/bin/bash

# should be run from the root of the project

NC=$(nproc --all)
DIR=$(mktemp -d /tmp/hfb3.XXXXXX)
INPUT=misc/tests/meta/input.hfb3

echo "===== Running in $DIR ====="

(make clean; make -j$NC SKILL=SKILL0) >/dev/null; bin/hfb3 $INPUT > $DIR/hfb3_SKILL0.out
echo SKILL0: $(cat $DIR/hfb3_SKILL0.out | grep -v SKILL | grep -v length | md5sum)

(make clean; make -j$NC SKILL=SKILL1) >/dev/null; bin/hfb3 $INPUT > $DIR/hfb3_SKILL1.out
echo SKILL1: $(cat $DIR/hfb3_SKILL1.out | grep -v SKILL | grep -v length | md5sum)

(make clean; make -j$NC SKILL=SKILL2) >/dev/null; bin/hfb3 $INPUT > $DIR/hfb3_SKILL2.out
echo SKILL2: $(cat $DIR/hfb3_SKILL2.out | grep -v SKILL | grep -v length | md5sum)

(make clean; make -j$NC SKILL=SKILL3) >/dev/null; bin/hfb3 $INPUT > $DIR/hfb3_SKILL3.out
echo SKILL3: $(cat $DIR/hfb3_SKILL3.out | grep -v SKILL | grep -v length | md5sum)

(make clean; make -j$NC SKILL=SKILL4) >/dev/null; bin/hfb3 $INPUT > $DIR/hfb3_SKILL4.out
echo SKILL4: $(cat $DIR/hfb3_SKILL4.out | grep -v SKILL | grep -v length | md5sum)

