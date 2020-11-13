#!/bin/csh
#
#PBS -S /bin/csh
#PBS -N test
#PBS -A NMMM0015
#PBS -l walltime=1:00
#PBS -q share
#PBS -o ./output_file
#PBS -j oe 
#PBS -k eod 
#PBS -l select=1:ncpus=6:mpiprocs=6
#PBS -m n    
#PBS -M schwartz@ucar.edu
#PBS -V 
#
#####PBS -m abe  # Main on abort, beginning, and ending (abe)
#####PBS -m n    # Don't mail
######PBS -M schwartz@ucar.edu
##########PBS -J 0-1000:10   (see http://www.arc.ox.ac.uk/content/pbs for "-J" option)
########PBS -l select=2:ncpus=36:mpiprocs=36:mem=110GB
########    set id = $PBS_ARRAYID     # would theoretically work with "-t"
#set id = $PBS_ARRAY_INDEX  # Works with -J

#
#BSUB -n 128
#BSUB -J hybrid
#BSUB -o output.driver
#BSUB -e output.driver
#BSUB -q regular
#BSUB -P NMMM0001
#BSUB -W 15
#BSUB -R "span[ptile=16]"

#touch -f ./hi_${id}

#touch -f test_done

cd $PWD

#setenv MPI_SHEPHERD true
mpiexec_mpt ./concatenate_netcdf_files.x
