#!/bin/bash

# Run n models for 100 time units, saving the last log line of each to
# file out to check energy conservation.

n=1
if (( $# > 0 )); then n=$1 ; fi
out=out
if (( $# > 1 )); then out=$2 ; fi

/bin/rm -f $out
touch $out
i=0
while (( $i < $n )) ; do
    makeplummer -n 5 -i -C | smallN -t 100 -d 1 | grep %% | tail -n 1 >> $out
    i=$((i+1))
done
