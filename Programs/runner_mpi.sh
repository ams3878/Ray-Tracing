#!/bin/bash
#

# This is an example bash script that is used to submit a job
# to the cluster.
#
# Typcially, the # represents a comment.However, #SBATCH is
# interpreted by SLURM to give an option from above. As you 
# will see in the following lines, it is very useful to provide
# that information here, rather than the command line.

# Name of the job - You MUST use a unique name for the job

# Standard out and Standard Error output files



# SBATCH -J c_sB_9
# SBATCH -o c_sB_9.out
# SBATCH -e c_sB_9.err

#SBATCH -J c_dB_T_T2
#SBATCH -o c_dB_T2.out
#SBATCH -e c_dB_T2.err

# In order for this to send emails, you will need to remove the
# space between # and SBATCH for the following 2 commands.
# Specify the recipient of the email
# SBATCH --mail-user=abc1234@rit.edu

# Notify on state change: BEGIN, END, FAIL or ALL
# SBATCH --mail-type=ALL

# Multiple options can be used on the same line as shown below.
# Here, we set the partition, number of cores to use, and the
# number of nodes to spread the jobs over.
#SBATCH -p kgcoe-mps -n 64 
#SBATCH --mem-per-cpu=4000M

#
# Your job script goes below this line.
#

# If the job that you are submitting is not sequential,
# then you MUST provide this line...it tells the node(s)
# that you want to use this implementation of MPI. If you
# omit this line, your results will indicate failure.
module load openmpi-1.10-x86_64

# Place your mpirun command here
# Notice that you have to provide the number of processes that
# are needed. This number needs to match the number of cores
# indicated by the -n option. If these do not, your results will
# not be valid or you may have wasted resources that others could
# have used.

# The following commands can be used as a starting point for 
# timing. Each example below needs to be modified to fit your parameters
# accordingly. You will also need to modify the config to point to
# box.xml to render the complex scene.
# **********************************************************************
# MAKE SURE THAT YOU ONLY HAVE ONE OF THESE UNCOMMENTED AT A TIME!
# **********************************************************************
# Sequential
# mpirun -np $SLURM_NPROCS raytrace_mpi -h 100 -w 100 -c configs/twhitted.xml -p none 
# Static Strips
# mpirun -np $SLURM_NPROCS raytrace_mpi -h 5000 -w 5000 -c configs/box.xml -p static_strips_horizontal

# Static Cycles
# mpirun -np $SLURM_NPROCS raytrace_mpi -h 5000 -w 5000 -c configs/box.xml -p static_cycles_vertical -cs 20

# Static Blocks
# mpirun -np $SLURM_NPROCS raytrace_mpi -h 5000 -w 5000 -c configs/box.xml -p static_blocks

# Dynamic
mpirun -np $SLURM_NPROCS raytrace_mpi -h 7000 -w 7000 -c configs/box.xml -p dynamic -bh 7 -bw 13