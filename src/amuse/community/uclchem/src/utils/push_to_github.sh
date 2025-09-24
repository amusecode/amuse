#!/usr/bin/bash

#not strictly necessary due to .gitignore but I like to clean up
cd src/fortran_src
make clean
cd ../..

#commit and push to github
git add . -A
#don't forget to update change log
git commit -F change.log
git push
