#!/bin/bash

#SBATCH --ntasks=1 
#SBATCH --time=12:00:00
#SBATCH --error=/sdf/home/c/cmantill/ldmx/pfreco/ldmx-sw/ldmx-fnal/batch/electron_inclusive/logs/slurm-%j.err
#SBATCH --output=/sdf/home/c/cmantill/ldmx/pfreco/ldmx-sw/ldmx-fnal/batch/electron_inclusive/logs/slurm-%j.out

source "$LDMX_BASE/ldmx-sw/scripts/ldmx-env.sh"
echo ${SLURM_ARRAY_JOB_ID}
echo ${SLURM_ARRAY_TASK_ID}

ldmx fire /sdf/home/c/cmantill/ldmx/pfreco/ldmx-sw/ldmx-fnal/batch/electron_inclusive/inclusive.py 50000 /sdf/home/c/cmantill/ldmx/pfreco/ldmx-sw/ldmx-fnal/batch/electron_inclusive/data/ ${SLURM_ARRAY_TASK_ID}
