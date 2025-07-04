#!/bin/bash
source /usr/share/Modules/init/bash

module load Programming_Languages/anaconda/3.11
conda activate /sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_clp
# Uncomment and adapt the following line
# if your kernel does not start
unset PYTHONPATH
export PYTHONPATH=/sps/lsst/groups/clusters/cl_pipeline_project/conda_envs/firecrown_clp:$PYTHONPATH
exec python -m ipykernel_launcher "$@" > debug.log 2>&1