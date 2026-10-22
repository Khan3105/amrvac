#!/bin/bash
# Run several WR_1D stages in sequence on the current node. Each stage after the first
# restarts from the last snapshot of the previous stage.
# Usage: bash output_revised/run_chain.sh <name1> [<name2> ...]
# (names = par files in output_revised/pars/ without .par = output folder names)
cd /home/r1082832/Codes/amrvac/tests/rhd/WR_1D || exit 1
export AMRVAC_DIR=/home/r1082832/Codes/amrvac       # needed at run time to find the OPAL tables
export OMP_NUM_THREADS=1
prev=""
for name in "$@"; do
  mkdir -p output_revised/$name
  out=output_revised/stdout/$name.out
  restart=""
  [ -n "$prev" ] && restart="-if $(ls output_revised/$prev/WR_1D_*.dat | tail -1)"
  echo "$(date) on $(hostname): starting $name $restart" | tee -a $out
  mpirun -np 8 --bind-to none ./output_revised/amrvac -i output_revised/pars/$name.par $restart 2>&1 | tee -a $out
  code=${PIPESTATUS[0]}
  echo "$(date): $name finished, exit code $code" | tee -a $out
  [ $code -ne 0 ] && break
  prev=$name
done
exec bash   # keep the screen window open so the end of the output stays visible
