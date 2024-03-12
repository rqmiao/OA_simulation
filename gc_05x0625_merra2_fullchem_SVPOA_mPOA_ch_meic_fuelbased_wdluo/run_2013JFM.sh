#! /bin/tcsh -f

# Set all PBS queue options, so no need to specify upon submit
# (retain hash-tag in front of all PBS command)
#PBS -N fuelbased_2013JFM
#PBS -q short
#PBS -e localhost:$PBS_O_WORKDIR/$PBS_JOBNAME.e$PBS_JOBID
#PBS -o localhost:$PBS_O_WORKDIR/$PBS_JOBNAME.o$PBS_JOBID
#PBS -l nodes=node03:ppn=40
#PBS -m ae
#PBS -M rqmiao@pku.edu.cn

# print some job info to the output file
cd $PBS_O_WORKDIR
echo `date`
echo Working directory is $PBS_O_WORKDIR
echo Running on host `hostname`
echo JobID is $PBS_JOBID

# set stacksize memory to unlimited
#limit stacksize unlimited

# change the stack size
#setenv KMP_STACKSIZE 100000000

# set correct ifort libraries
#source /share/apps/intel/bin/ifortvars.csh intel64

# Set files
set inp = "input.geos.2013JFM"
set hemco = "HEMCO_Config_run.rc"
set history = "HISTORY_run.rc"

# Copy input file
if ( -f $inp ) then
   cp -f $inp ./input.geos
   cp -f $hemco ./HEMCO_Config.rc
   cp -f $history ./HISTORY.rc
else
   echo "Could not find input file!"
   exit(1)
endif

# Remove previous log and run the code
rm -f log_$PBS_JOBNAME
time ./gcclassic > log_$PBS_JOBNAME

# Exit normally
exit(0)
