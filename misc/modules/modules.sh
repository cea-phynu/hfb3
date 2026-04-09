# module config for inti-amd.ocre.cea.fr

module --force purge || true

module load gnu/11 mpi/openmpi/4 armadillo/15.2.3 python3/3.12 lapack/netlib boost

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
