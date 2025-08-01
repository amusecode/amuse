#!/bin/bash

work_dir=multiples_tmp

N=1000
salpeter="-s"

infile='t=0000.00'
t=20
d=1.0
D=5.0
savemode='-S'

#rm -rf $work_dir
if [ ! -e $work_dir ]; then
    mkdir $work_dir
fi

cd $work_dir
$AMUSE_DIR/amuse.sh ../initialize_system.py -n $N $salpeter
$AMUSE_DIR/amuse.sh ../run_ph4.py -i $infile -d $d -D $D -t $t $savemode \
		> t=0.log
