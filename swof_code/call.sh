#!/bin/bash
### Job Name
#PBS -N mpi_job
### Project code
#PBS -A UDKE0018
#PBS -l walltime=01:00:00
#PBS -q queue_name
### Merge output and error files
#PBS -j oe
#PBS -k eod
### Select 2 nodes with 36 CPUs each for a total of 72 MPI processes
#PBS -l select=2:ncpus=36:mpiprocs=36
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M quatratavia@gmail.com

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run the executable
python call_james_ncar.py