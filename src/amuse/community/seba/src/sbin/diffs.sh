#!/bin/bash

for file in `find include -name \*.\[h\]` ; do
    echo comparing $file...
    diff -bwi $file $STARLAB/$file
done

for file in `find {std,node} -name \*.\[Cc\]` ; do
    echo comparing $file...
    diff -bwi $file $STARLAB/src/$file
done

for file in `find sstar -name \*.\[Cc\]` ; do
    echo comparing $file...
    diff -bwi $file $STARLAB/src/star/$file
done
