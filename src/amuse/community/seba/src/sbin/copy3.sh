#!/bin/bash

for file in `find sstar -name \*.\[Cc\]` ; do
    echo comparing $file...
    diff -bwi $file $STARLAB/src/star/$file
done
