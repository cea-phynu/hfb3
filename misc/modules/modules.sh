# module config for inti-amd.ocre.cea.fr

module --force purge || true

module load ccc
module load dfldatadir
module load mkl
module load armadillo
module load python3
module load swig
module load gnu
module load boost

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
