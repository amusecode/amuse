#!/bin/bash

RAMSES=`dirname $0`/../../..
echo "Running from ramses directory ${RAMSES}"

mkdir test
cd test
echo "Storing results in" `pwd`


# This takes several minutes.
cp $RAMSES/aton/sims/testing/stromgren.nml .
cp $RAMSES/aton/sims/testing/expected_profile.txt .
mpirun -n 1 $RAMSES/bin/ramses3d stromgren.nml | tee log

# Get some interesting data.
$RAMSES/utils/f90/amr2cell -inp output_00007/ -out output_00007/cells.txt
python $RAMSES/aton/utils/test5spherical.py < output_00007/cells.txt > profile.txt

cd ..

echo ""
echo ""
echo "********************************************************************************"
echo "Now use gnuplot to compare test/profile.txt to sims/testing/expected_profile.txt"
echo "The columns are r, density, xneutral, pressure, temperature, mach, xion"
echo "e.g. to compare the temperature in gnuplot:"
echo "set log y"
echo "plot 'test/profile.txt' using 1:5, 'test/expected_profile.txt' using 1:5"

