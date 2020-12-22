#!/bin/sh

## how many MPI tasks in total
#SBATCH --ntasks=1

## how many CPUs for each task
#SBATCH --cpus-per-task=6

## how many nodes will use in total
#SBATCH --nodes=1

## how to distribute the tasks among nodes
#SBATCH --ntasks-per-node=1


#SBATCH -t 3:00:00
#SBATCH -p bigmem
#SBATCH --mem=512G
#SBATCH -J QOSF
##SBATCH -e slurm-%x-%j.err
#SBATCH --output=slurm-%x-%j.out

module load python/3.6.6 gcc/8.3
cd ~/scratch/QOSF-FeMoco2020/
python3 ./src/gen_hamiltonian.py

