#!/bin/sh

echo $1
readelf -d  $1 | grep RPATH
