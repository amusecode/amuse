#!/bin/bash

for file in `find {std,node} -name \*.\[Cc\]` ; do
    cp -p $STARLAB/src/$file $file
done
