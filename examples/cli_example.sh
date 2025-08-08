#!/bin/bash

# This example shows how to chain some HFB calculations.
# It should be launched from the root of the project.

# generate `bin/hfb3` if needed
make -j4

# use only one CPU core
export OMP_NUM_THREADS=1

# do a first calculation for 16O with <q20>=8fm^2 and <q30>=0fm^3

echo '#######################'
echo '## first calculation ##'
echo '#######################'

bin/hfb3 examples/16O_deformed.hfb3

# calculate some observables of the resulting state

echo '#################'
echo '## observables ##'
echo '#################'

bin/hfb3 16O_deformed.msg.gz '{action:obs}'

# do a second calculation for 16O with <q20>=10fm^2, no constraint on <q30>, no basis optimization, starting from the previous result

echo '########################'
echo '## second calculation ##'
echo '########################'

bin/hfb3 16O_deformed.msg.gz '{constraints/q20t:10.0,constraints/q30t:,action/basisOptimization:false,action/jobName:16O_q20t_10}'

# calculate some observables of the resulting state

echo '#################'
echo '## observables ##'
echo '#################'

bin/hfb3 16O_q20t_10.msg.gz '{action:obs}'


