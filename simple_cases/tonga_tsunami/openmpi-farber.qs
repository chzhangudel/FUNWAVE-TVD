#
# Template:  OpenMPI, Low-Bandwidth (Gigabit Ethernet) Variant
# Revision:  $Id: openmpi-gige.qs 577 2015-12-21 14:39:43Z frey $
#
# Usage:
# 1. Modify "NPROC" in the -pe line to reflect the number
#    of processors desired.
# 2. Modify the value of "MY_EXE" to be your MPI program and 
#    "MY_EXE_ARGS" to be the array of arguments to be passed to
#    it.
# 3. Uncomment the WANT_CPU_AFFINITY line if you want Open MPI to
#    bind workers to processor cores.  Can increase your program's
#    efficiency.
# 4. Uncomment the SHOW_MPI_DEBUGGING line if you want very verbose
#    output written to the Grid Engine output file by OpenMPI.
# 5. If you use exclusive=1, please be aware that NPROC will be
#    rounded up to a multiple of 20.  In this case, set the
#    WANT_NPROC variable to the actual core count you want.  The
#    script will "load balance" that core count across the N nodes
#    the job receives.
# 6. Jobs default to using 1 GB of system memory per slot.  If you
#    need more than that, set the m_mem_free complex.
# 
#$ -pe mpi 20
#
# Change the following to #$ and set the amount of memory you need
# per-slot if you're getting out-of-memory errors using the
# default:
# -l m_mem_free=2G
#
# Change the following to #$ if you want exclusive node access
# (see 5. above):
#$ -l exclusive=1
#$ -l standby=0
#
# If you want an email message to be sent to you when your job ultimately
# finishes, edit the -M line to have your email address and change the
# next two lines to start with #$ instead of just #
#$ -m eas
#$ -M chzhang@udel.edu
#

#
# Setup the environment; choose the OpenMPI version that's
# right for you:
#
vpkg_devrequire openmpi/1.8.2-intel64

#
# The MPI program to execute:
#
MY_EXE="/home/1670/src/src_FUNWAVEM/funwave-METEO--mpif90-parallel-double"

#
# Arguments to the MPI program being executed.  Remember to use quotes
# around any arguments with whitespace or special characters, e.g.
#
#   MY_EXE_ARGS=("this is arg1" arg2 arg3)
MY_EXE_ARGS=()

#
# By default the slot count granted by Grid Engine will be
# used, one MPI worker per slot.  Set this variable if you
# want to use fewer cores than Grid Engine granted you (e.g.
# when using exclusive=1):
#
#WANT_NPROC=0

#
# Ask Open MPI to do processor binding?
#
#WANT_CPU_AFFINITY=YES

#
# Uncomment to enable lots of additional information as OpenMPI executes
# your job:
#
#SHOW_MPI_DEBUGGING=YES

##
## You should NOT need to change anything after this comment.
##
OPENMPI_FLAGS="--display-map --mca btl ^openib --mca oob_tcp_if_exclude lo,ib0 --mca btl_tcp_if_exclude lo,ib0"
if [ "${WANT_CPU_AFFINITY:-NO}" = "YES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --bind-to core"
fi
if [ "${WANT_NPROC:-0}" -gt 0 ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --np ${WANT_NPROC} --map-by node"
fi
if [ "${SHOW_MPI_DEBUGGING:-NO}" = "YES" ]; then
  OPENMPI_FLAGS="${OPENMPI_FLAGS} --debug-devel --debug-daemons --display-devel-map --display-devel-allocation --mca mca_verbose 1 --mca coll_base_verbose 1 --mca ras_base_verbose 1 --mca ras_gridengine_debug 1 --mca ras_gridengine_verbose 1 --mca btl_base_verbose 1 --mca mtl_base_verbose 1 --mca plm_base_verbose 1 --mca pls_rsh_debug 1"
  if [ "${WANT_CPU_AFFINITY:-NO}" = "YES" ]; then
    OPENMPI_FLAGS="${OPENMPI_FLAGS} --report-bindings"
  fi
fi

echo "GridEngine parameters:"
echo "  mpirun        = "`which mpirun`
echo "  nhosts        = $NHOSTS"
echo "  nproc         = $NSLOTS"
echo "  executable    = $MY_EXE"
echo "  MPI flags     = $OPENMPI_FLAGS"
echo "-- begin OPENMPI run --"
mpirun ${OPENMPI_FLAGS} "$MY_EXE" "${MY_EXE_ARGS[@]}"
rc=$?
echo "-- end OPENMPI run --"

exit $rc

