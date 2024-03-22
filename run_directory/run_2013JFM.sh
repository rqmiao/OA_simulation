#! /bin/tcsh -f

# Set all PBS queue options, so no need to specify upon submit
# (retain hash-tag in front of all PBS command)
#PBS -N OA_2013JFM
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

# Remove previous log and run the code
rm -f log_$PBS_JOBNAME
time ./gcclassic > log_$PBS_JOBNAME

# Exit normally
exit(0)
