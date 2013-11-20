#! /bin/bash

bin=$HOME'/bin'

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
lbox=`grep unit_l output_$count/info_* | cut -d = -f 2`
h0=`grep H0 output_$count/info_* | cut -d = -f 2`
echo 'aexp='$aexp
echo 'lbox='$lbox
echo 'h0='$h0
redshift=`echo $aexp|awk '{OFMT="%.4f"; print 1/$1-1}'`
hh=`echo $h0|awk '{OFMT="%.4f"; print $1/100.}'`
ll=`echo $lbox|awk '{OFMT="%.4f"; print $1/3.08E+21}'`
echo 'z='$redshift
echo 'h='$hh
echo 'l='$ll

source params_$count.sh

r200=`echo $rvir $ll|awk '{print (($1)/($2))}'`
rgal=`echo $r200|awk '{print $1/10}'`
hgal=`echo $rgal|awk '{print $1/10}'`

echo 'rvir='$rvir
echo 'r200c='$r200
echo 'rgal='$rgal
echo 'hgal='$hgal

echo $xc
echo $yc
echo $zc

xmi=`echo $xc $imsize|awk '{print $1-$2}'`
xma=`echo $xc $imsize|awk '{print $1+$2}'`
ymi=`echo $yc $imsize|awk '{print $1-$2}'`
yma=`echo $yc $imsize|awk '{print $1+$2}'`
zmi=`echo $zc $imsize|awk '{print $1-$2}'`
zma=`echo $zc $imsize|awk '{print $1+$2}'`

part2prof -inp output_$count -out virial_$count.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -rma $r200 -nra 400

part2prof -inp output_$count -out gal_$count.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -rma $rgal -nra 400

if [ -z "$2" ]; then
    exit
fi

if [ $2 eq 0 ]; then
    exit
fi

amr2prof  -inp output_$count -out virial_$count.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -rma $r200 -nra 400

paste virial_$count.prof.dark virial_$count.prof.star virial_$count.prof.gas > virial_$count.prof.tot

part2prof -inp output_$count -out gal_$count.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -rma $rgal -nra 400

amr2prof  -inp output_$count -out gal_$count.prof -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -rma $rgal -nra 400

paste gal_$count.prof.dark gal_$count.prof.star gal_$count.prof.gas > gal_$count.prof.tot

amr2cylprof -inp output_$count -out gal_$count.cyl -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -jx $jxc -jy $jyc -jz $jzc -rma $rgal -nra 200 -hma $hgal

part2cylprof -inp output_$count -out gal_$count.cyl -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -jx $jxc -jy $jyc -jz $jzc -rma $rgal -nra 200 -hma $hgal -cir  gal_$count.prof.tot

part2cylprof -inp output_$count -out virial_$count.cyl -xce $xc -yce $yc -zce $zc -uce $uc -vce $vc -wce $wc -jx $jxc -jy $jyc -jz $jzc -rma $r200 -nra 400 -hma $r200 -cir  virial_$count.prof.tot

