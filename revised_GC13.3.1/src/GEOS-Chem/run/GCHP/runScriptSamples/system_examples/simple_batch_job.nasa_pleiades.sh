#!/bin/bash
#
# DESCRIPTION      Sample GCHP batch job script for Pleiades
# SCHEDULER        PBS
# DATE             2021-02-18
# SIMULATION       C90, 2019-01-01 to 2019-01-02
# RESOURCES        2 nodes (Broadwell; model=bro), 24 processes per node
# SEE ALSO         qsub man pages, http://gchp.rtfd.io/
#
#PBS -S /bin/bash
#PBS -l select=2:ncpus=24:mpiprocs=24:model=bro
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -W group_list=xXXXX
#PBS -m e

# Load modules (so software dependencies are available)
module purge
module use /nobackup/lbindle/modulefiles
module load gchp-env/2021.06-gnu

# Misc. configuration
module list     # print loaded modules
set -e          # if a subsequent command fails, treat it as fatal (don't continue)
set -x          # for remainder of script, echo commands to the job's log file

# Unlimit resources (to prevent OS killing GCHP due to resource usage)
ulimit -c 0                  # coredumpsize
ulimit -l unlimited          # memorylocked
ulimit -u 50000              # maxproc
ulimit -v unlimited          # vmemoryuse
ulimit -s unlimited          # stacksize

# cd to working directory (by default, PBS jobs land in $HOME)
cd $PBS_O_WORKDIR

# Simple GCHP launch
rm -f cap_restart                   # delete restart start time file if present
./runConfig.sh                      # update configuration files
mpiexec -n 48 ./gchp > runlog.txt   # launch 48 GCHP processes

