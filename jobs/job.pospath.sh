#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --output=slurm.out
#SBATCH --job-name=%%jobname%%
#SBATCH --nodes=1
#SBATCH --tasks-per-node=40
#SBATCH --mem=187000

#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
##SBATCH --constraint=IB

################################################################################
# Here is a example for UB CCR.
# Modify as your need.


cd $SLURM_SUBMIT_DIR
date
module load intel mkl intel-mpi
module list
ulimit -s  unlimited

echo "Determining node information ..."
SLURM_NODEFILE=nodes.$SLURM_JOB_ID
srun hostname -s > $SLURM_NODEFILE

export NODELIST=tmp.$$
echo "group main" > $NODELIST
NODES=`cat $SLURM_NODEFILE`
for node in $NODES ; do
   echo "host "$node >> $NODELIST
done

echo "NODELIST = "$NODELIST
echo "nodes = $SLURM_JOB_NUM_NODES"
PROCS_ON_NODE=$((SLURM_NPROCS/SLURM_JOB_NUM_NODES))
echo "tasks_per_node = $PROCS_ON_NODE"
echo "total_procs = $SLURM_NPROCS"
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

################################################################################

MPIRUN="srun --nodes=$SLURM_JOB_NUM_NODES --ntasks-per-node=$PROCS_ON_NODE --cpus-per-task=1"

VASP_DIR="/projects/academic/pzhang3/zhao_code/vasp.6.1.0"
VASP="$VASP_DIR/bin/vasp_std"
VASPG="$VASP_DIR/bin/vasp_gam"
VASPZ="$VASP_DIR/bin/vasp_std_z"

################################################################################

ln -s ../POTCAR
ln -s ../INCAR
ln -s ../KPOINTS
$MPIRUN $VASPG > vasp.out

################################################################################

rm tmp.$$
rm nodes.$SLURM_JOB_ID
echo "All Done!"
date
