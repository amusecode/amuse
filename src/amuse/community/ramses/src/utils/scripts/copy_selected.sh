#!/bin/bash
mkdir -p $2
while read line; do 
    cp -r $line $2
done < $1
