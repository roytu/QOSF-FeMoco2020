#!/bin/sh

## how many MPI tasks in total
#SBATCH --ntasks=1

## how many CPUs for each task
#SBATCH --cpus-per-task=8

## how many nodes will use in total
#SBATCH --nodes=1

## how to distribute the tasks among nodes
#SBATCH --ntasks-per-node=1


#SBATCH -t 12:00:00
#SBATCH -p bigmem
#SBATCH --mem=512G
#SBATCH -J QOSF
##SBATCH -e slurm-%x-%j.err
#SBATCH --output=slurm-%x-%j.out

module load python/3.9.0 gcc/8.3 cmake/3.8.0
cd ~/scratch/QOSF-FeMoco2020/
source ~/scratch/qiskitenv/bin/activate
python3 ./src/gen_hamiltonian_qiskit_test.py

