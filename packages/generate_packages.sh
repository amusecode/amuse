#!/bin/bash

mkdir -p dist

for p in amuse*
do
  cd $p || exit $?
  rm -rf dist || exit $?
  python setup.py sdist || exit $?
  cp dist/*.tar.gz ../dist/ || exit $?
  cd .. 
done
