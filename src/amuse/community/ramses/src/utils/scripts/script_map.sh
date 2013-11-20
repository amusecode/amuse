#! /bin/bash

type='cdm'
dmin=1e-10
dmax=1e-6

bin=$HOME'/bin'

a=$1
if [ $a -ge 100 ]; then
    count='00'$a
else
    count='000'$a
fi
echo output_$count
i=$(($a-10))

r0=0.473753
r1=0.000297422
r2=-1.06913e-06
r3=1.46028e-09
xc=`echo $i $r0 $r1 $r2 $r3 |awk '{print $2+$3*$1+$4*$1*$1+$5*$1*$1*$1}'`

r0=0.547496
r1=-0.000548908
r2=1.76895e-06
r3=-2.56450e-09
yc=`echo $i $r0 $r1 $r2 $r3 |awk '{print $2+$3*$1+$4*$1*$1+$5*$1*$1*$1}'`

r0=0.505470
r1=-0.000101028
r2=4.91462e-07
r3=-6.61523e-10
zc=`echo $i $r0 $r1 $r2 $r3 |awk '{print $2+$3*$1+$4*$1*$1+$5*$1*$1*$1}'`

aexp=`grep aexp output_$count/info_* | cut -d = -f 2`
imsize=`echo $aexp|awk '{print 0.0025/$1}'` 
echo 'aexp='$aexp
xmi=`echo $xc $imsize|awk '{print $1-$2}'`
xma=`echo $xc $imsize|awk '{print $1+$2}'`
ymi=`echo $yc $imsize|awk '{print $1-$2}'`
yma=`echo $yc $imsize|awk '{print $1+$2}'`
zmi=`echo $zc $imsize|awk '{print $1-$2}'`
zma=`echo $zc $imsize|awk '{print $1+$2}'`

$bin/part2map -inp output_$count -out cdm_map_$count.dat -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 256 -ny 256 -dir y
# -str true -age true
$bin/map2img.py cdm_map_$count.dat -l -c jet

#$bin/amr2map -inp output_$count -out dens_map_$count.dat -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -lma 16 -nx 256 -ny 256 -dir z -typ 1 
#$bin/map2img.py dens_map_$count.dat -l -c jet

