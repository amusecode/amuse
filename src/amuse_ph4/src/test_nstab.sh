#!/bin/bash

# Check that the Fortran and C versions of nstab actually return the
# same results!  Run n random tests.

n=100
if (( $# > 0 )); then n=$1 ; fi

make -f Makefile.ph4 test_nstab

i=0
bad=0
while (( $i < $n )); do
    i=$((i+1))
    s=`./random_nstab`
    f=`ftest_nstab $s`
    c=`ctest_nstab $s`
    if [ $f != $c ]; then
	bad=$(($bad+1))
	echo $s `ftest_nstab $s` `ctest_nstab $s`
    fi
done

echo $bad disagreement\(s\) found
