#!/bin/sh

# Compile BIOS

rm ./bios
printf "Compiling BIOS...\n\n"
make -s

# Optimization Data File (.opt):

FILE=./Args

# Hostfile:

HOSTS=./mpi_hosts

printf "BIOS will run on the following nodes:\n";
cat $HOSTS

printf "\nOptimization model file to be processed: $FILE.opt \n\n"

# Number of MPI processes

NUMPROC=5

printf "The number of MPI Processes is: $NUMPROC \n\n";

# BIOS folder:

BIOSFOLDER=./Desktop/BIOSv2Public/BIOSv2Public/

# Copy folder:

COPYFOLDER=/tmp/BIOS

# Copy BIOS to each node.

printf "Copying BIOS to each node...\n\n";

pdsh -a rm -rf $COPYFOLDER

pdsh -a cp -r $BIOSFOLDER $COPYFOLDER

# Run BIOS in each node.

printf "Running optimization process...\n\n"

/usr/local/bin/mpirun -wd $COPYFOLDER -hostfile $HOSTS -np $NUMPROC $COPYFOLDER/bios $FILE

# Copy results back to original folder

#ssh n001 scp $COPYFOLDER/$FILE.out iuri@clusterlmcv:$BIOSFOLDER/$FILE-n001.out
#ssh n002 scp $COPYFOLDER/$FILE.out iuri@clusterlmcv:$BIOSFOLDER/$FILE-n002.out
#ssh n003 scp $COPYFOLDER/$FILE.out iuri@clusterlmcv:$BIOSFOLDER/$FILE-n003.out
#ssh n004 scp $COPYFOLDER/$FILE.out iuri@clusterlmcv:$BIOSFOLDER/$FILE-n004.out
#ssh n005 scp $COPYFOLDER/$FILE.out iuri@clusterlmcv:$BIOSFOLDER/$FILE-n005.out

printf "Optimization process ended.\n"
