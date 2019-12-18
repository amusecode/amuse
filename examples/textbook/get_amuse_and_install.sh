#! /bin/bash

GITHASH=0b7d0c5aa8dcaa2fc564aba62dc9a7f686c29382

temp=amuse-temp-$$
mkdir $temp; cd $temp

git clone https://github.com/amusecode/amuse.git
cd amuse; git reset --hard $GITHASH; cd ../..

export AMUSE_DIR="$PWD/amuse-$GITHASH"
mv $temp/amuse $AMUSE_DIR
rm -rf $temp

cd $AMUSE_DIR
./configure && make framework
