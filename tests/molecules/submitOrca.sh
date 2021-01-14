#!/bin/bash
#SBATCH --time=200:00:00
#SBATCH -n 32

#load modules
module load orca/4.2.1 mpi/openmpi_4.0.4_gcc gcc/8.3

#Execute
cd ~/QOSF-FeMoco2020/molecules
export RSH_COMMAND="/usr/bin/ssh -x"
export ORCA_DIR="/gpfs/runtime/opt/orca/4.2.1/bin"
$ORCAMP/orca ~/QOSF-FeMoco2020/molecules/$1 > ~/QOSF-FeMoco2020/molecules/${1%.inp}.out

