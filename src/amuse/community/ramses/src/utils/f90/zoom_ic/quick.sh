#! /bin/bash

#bin='/project/s201/martdav/ramses/utils/f90/zoom_ic/'

cd /scratch/rosa/martdav/ref_box/dmo/

a=$1
count='0000'$a
if [ $a -ge 10 ]; then
    count='000'$a
fi
if [ $a -ge 100 ]; then
    count='00'$a
fi
echo output_$count
aexp=`grep aexp output_$count/info_* | cut -d = -f 2`
echo 'aexp='$aexp
redshift=`echo $aexp|awk '{OFMT="%.2f"; print 1/$1-1}'`
echo 'z='$redshift

dir='y'

xc=$2
yc=$3
zc=$4
rvir=$5
number=$6

halocount='0000'$number
if [ $number -ge 10 ]; then
    halocount='000'$number
fi
if [ $number -ge 100 ]; then
    halocount='00'$number
fi
if [ $number -ge 1000 ]; then
    halocount='0'$number
fi

rad=`echo $rvir|awk '{print $1*2.5}'`

ic=`echo $xc|awk '{print $1*512}'`
jc=`echo $yc|awk '{print $1*512}'`
kc=`echo $zc|awk '{print $1*512}'`

xmi=`echo $xc $imsize|awk '{print $1-$2}'`
xma=`echo $xc $imsize|awk '{print $1+$2}'`
ymi=`echo $yc $imsize|awk '{print $1-$2}'`
yma=`echo $yc $imsize|awk '{print $1+$2}'`
zmi=`echo $zc $imsize|awk '{print $1-$2}'`
zma=`echo $zc $imsize|awk '{print $1+$2}'`

#$bin/geticref -xc $xc -yc $yc -zc $zc -rad $rad -inp output_$count -per true
#$bin/geticmask -inp output_00001 -smt 4 -rsm 12 -gfc ../../ic_files/ref_box_512 -per true

geticref -xc $xc -yc $yc -zc $zc -rad $rad -inp output_$count -per true
geticmask -inp output_00001 -smt 4 -rsm 12 -gfc ../../ic_files/ref_box_512 -per true

cd ../../ic_halos

mkdir halo_$halocount

cd halo_$halocount

mkdir ic_orig_512
mkdir ic_b200_512
mkdir ic_b200_256
mkdir ic_b200_128
mkdir ic_b200_64
mkdir ic_b100_256
mkdir ic_b100_128

cd ic_orig_512
ln -s ../../../ic_files/ref_box_512/ic_deltab 
ln -s ../../../ic_files/ref_box_512/ic_velcx
ln -s ../../../ic_files/ref_box_512/ic_velcy
ln -s ../../../ic_files/ref_box_512/ic_velcz
ln -s ../../../ic_files/ref_box_512/ic_velbx
ln -s ../../../ic_files/ref_box_512/ic_velby
ln -s ../../../ic_files/ref_box_512/ic_velbz
ln -s ../../../ref_box/dmo/ic_refmap
#ln -s ../../../ref_box/dmo/ic_pvar_00001
cd ..

cat > inp.graf <<EOF
$ic $jc $kc
EOF

#$bin/center_grafic ic_orig_512 ic_b200_512 < inp.graf

center_grafic ic_orig_512 ic_b200_512 < inp.graf
rm -rf ic_orig_512

#$bin/degrade_grafic ic_b200_512 ic_b200_256
#$bin/degrade_grafic ic_b200_256 ic_b200_128
#$bin/degrade_grafic ic_b200_128 ic_b200_64

degrade_grafic ic_b200_512 ic_b200_256
degrade_grafic ic_b200_256 ic_b200_128
degrade_grafic ic_b200_128 ic_b200_64

cat > inp.graf <<EOF
256 256 256
256 256 256
EOF

#$bin/extract_grafic ic_b200_512 ic_b100_256 < inp.graf
#$bin/degrade_grafic ic_b100_256 ic_b100_128

extract_grafic ic_b200_512 ic_b100_256 < inp.graf
degrade_grafic ic_b100_256 ic_b100_128

rm inp.graf

cat > tarb.sh <<EOF
tar -cf ic_b200_512.tar ic_b200_512
tar -cf ic_b200_256.tar ic_b200_256
tar -cf ic_b200_128.tar ic_b200_128
tar -cf ic_b200_64.tar ic_b200_64
tar -cf ic_b100_256.tar ic_b100_256
tar -cf ic_b100_128.tar ic_b100_128
EOF

sh tarb.sh

cd ../../run_halos/

mkdir halo_$halocount 

cd halo_$halocount 

cat > tara.sh <<EOF
ls -d output_* > list_o
for i in `cat list_o`; do tar -cf $i.tar $i; done
rm list_o
EOF

cat > dmo.nml <<EOF

&RUN_PARAMS
cosmo=.true.
pic=.true.
poisson=.true.
hydro=.false.
nrestart=0
nremap=1
nsubcycle=1,1,1,2,2,2,2 
ncontrol=1
/

&OUTPUT_PARAMS
delta_aout=0.01
aend=1.0
/

&INIT_PARAMS
filetype='grafic'
initfile(1)='/scratch/rosa/martdav/ic_halos/halo_$halocount/ic_b200_128'
initfile(2)='/scratch/rosa/martdav/ic_halos/halo_$halocount/ic_b200_256'
initfile(3)='/scratch/rosa/martdav/ic_halos/halo_$halocount/ic_b100_256'
/

&AMR_PARAMS
levelmin=7
levelmax=19
ngridmax=1000000
npartmax=4000000
nexpand=4,1,1,1,1,1
/

&POISSON_PARAMS
epsilon=1.d-3
/

&REFINE_PARAMS
m_refine=20*8.,
ivar_refine=0
mass_cut_refine=1e-8
/

EOF


cat > job_ramses.pbs <<EOF
#! /bin/bash
#PBS -l mppwidth=64
#PBS -l walltime=24:00:00
#PBS -V

set -ex

cd $SCRATCH/run_halos/halo_$halocount

export DATE=`date +%F_%Hh%M`

aprun -n 64 /project/s201/martdav/ramses/bin/ramses3d dmo.nml > run$DATE.log

exit
EOF

#qsub job_ramses.pbs
#sh tara.sh

cd /scratch/rosa/martdav/ref_box/dmo/
rm ic_ref*
