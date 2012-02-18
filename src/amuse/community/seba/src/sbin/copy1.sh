#!/bin/bash

for file in `find include -name \*.\[h\]` ; do
    cp -p $STARLAB/$file $file
done
