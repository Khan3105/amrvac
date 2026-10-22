#!/bin/bash
# Start one revised WR_1D run on the current node. Usage: bash output_revised/run.sh <run name>
# (the run name is the par file name in output_revised/pars/ without .par, and the output folder name)
cd /home/r1082832/Codes/amrvac/tests/rhd/WR_1D || exit 1
export AMRVAC_DIR=/home/r1082832/Codes/amrvac       # needed at run time to find the OPAL tables
export OMP_NUM_THREADS=1
name=$1
mkdir -p output_revised/$name
out=output_revised/stdout/$name.out
echo "$(date) on $(hostname): starting $name" | tee -a $out
mpirun -np 8 --bind-to none ./output_revised/amrvac -i output_revised/pars/$name.par 2>&1 | tee -a $out
echo "$(date): $name finished, exit code ${PIPESTATUS[0]}" | tee -a $out
exec bash   # keep the screen window open so the end of the output stays visible
