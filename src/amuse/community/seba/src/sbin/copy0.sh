#!/bin/bash

for file in `find include -name config.h` ; do
    cp -p $STARLAB/$file $file
done
