#!/bin/sh
### Set the job name (for your reference)
#PBS -N COL380_A3
### Set the project name, your department code by default
#PBS -P cse
### Request email when job begins and ends
#PBS -m bea
### Specify email address to use for notification.
### #PBS -M $USER@iitd.ac.in
#### FORMAT: -l select=n:ncpus=m
#PBS -l select=1:ncpus=2:ngpus=0:mem=4G
### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=00:01:00

#PBS -l software=c++
# After job starts, must goto working directory. 
# $PBS_O_WORKDIR is the directory from where the job is fired.
set "ARGS311=--taskid=1 --verbose=1 --startk=3 --endk=10 --inputpath=A3/test1/test-input-1.gra --headerpath=A3/test1/test-header-1.dat --outputpath=A3/test1/task1-output-1.txt"
set "ARGS312 = --taskid=1 --verbose=1 --startk=3 --endk=12 --inputpath=A3/test2/test-input-2.gra --headerpath=A3/test2/test-header-2.dat --outputpath=A3/test2/task1-output-2.txt"
set "ARGS313 = --taskid=1 --verbose=1 --startk=4 --endk=8 --inputpath=A3/test3/test-input-3.gra --headerpath=A3/test3/test-header-3.dat --outputpath=A3/test3/task1-output-3.txt"
set "ARGS314 = --taskid=1 --verbose=1 --startk=2 --endk=3 --inputpath=A3/test4/test-input-4.gra --headerpath=A3/test4/test-header-4.dat --outputpath=A3/test4/task1-output-4.txt"
# A3: TASK 2 Test Cases
set "ARGS321 = --taskid=2 --verbose=1 --endk=3 --inputpath=A3/test1/test-input-1.gra --headerpath=A3/test1/test-header-1.dat --outputpath=A3/test1/task2-output-1.txt --p=10"
set "ARGS322 = --taskid=2 --verbose=1 --startk=1 --endk=3 --inputpath=A3/test2/test-input-2.gra --headerpath=A3/test2/test-header-2.dat --outputpath=A3/test2/task2-output-2.txt --p=10"
set "ARGS323 = --taskid=2 --verbose=1 --startk=1 --endk=3 --inputpath=A3/test3/test-input-3.gra --headerpath=A3/test3/test-header-3.dat --outputpath=A3/test3/task2-output-3.txt --p=4"
set "ARGS324 = --taskid=2 --verbose=1 --startk=1 --endk=2 --inputpath=A3/test4/test-input-4.gra --headerpath=A3/test4/test-header-4.dat --outputpath=A3/test4/task2-output-4.txt --p=20"

module load compiler/cuda/9.2/compilervars
module load compiler/gcc/9.1.0
module load compiler/gcc/6.5/openmpi/4.0.2
module load compiler/gcc/9.1/openmpi/4.0.2
module load compiler/gcc/9.1/mpich/3.3.1

echo "==============================="
echo $PBS_JOBID
cat $PBS_NODEFILE
echo "==============================="
cd $PBS_O_WORKDIR
#job 
mpirun -n {n*m} ./a3 --taskid=1 --verbose=1 --startk=3 --endk=10 --inputpath=A3/test1/test-input-1.gra --headerpath=A3/test1/test-header-1.dat --outputpath=A3/test1/task1-output-1.txt
#NOTE
# The job line is an example : users need to change it to suit their applications
# The PBS select statement picks n nodes each having m free processors
# OpenMPI needs more options such as $PBS_NODEFILE
