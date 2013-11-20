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
echo 'aexp='$aexp
redshift=`echo $aexp|awk '{OFMT="%.2f"; print 1/$1-1}'`
echo 'z='$redshift

dir='y'
xc=0.499055
yc=0.500313
zc=0.500152
imsize=0.01
lma=17

xmi=`echo $xc $imsize|awk '{print $1-$2}'`
xma=`echo $xc $imsize|awk '{print $1+$2}'`
ymi=`echo $yc $imsize|awk '{print $1-$2}'`
yma=`echo $yc $imsize|awk '{print $1+$2}'`
zmi=`echo $zc $imsize|awk '{print $1-$2}'`
zma=`echo $zc $imsize|awk '{print $1+$2}'`


part2map -inp output_$count -out cdm_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 512 -ny 512 -dir ${dir} -den hop/hop$count.den #-fil ascii 

sunset -inp output_$count -out i_band_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 512 -ny 512 -dir ${dir} -bnd i_prime 

#./part2map -inp output_$count -out cdm.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 1024 -ny 1024 -dir ${dir}

#part2map -inp output_$count -out star_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -nx 512 -ny 512 -dir ${dir} -str true -age true -fil ascii

amr2map -inp output_$count -out dens_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -lma ${lma} -nx 512 -ny 512 -typ 1 -dir ${dir} #-fil ascii
amr2map -inp output_$count -out temp_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -lma $lma -nx 512 -ny 512 -typ 5 -dir ${dir} #-fil ascii
amr2map -inp output_$count -out metal_${count}_dir${dir}.map -xmi $xmi -xma $xma -ymi $ymi -yma $yma -zmi $zmi -zma $zma -lma $lma -nx 512 -ny 512 -typ 6 -dir ${dir} #-fil ascii

map2img.py cdm_${count}_dir${dir}.map -o cdm_${count}_dir${dir}.png -l -c gray -m 1e1 -M 1e7
map2img.py i_band_${count}_dir${dir}.map -o i_band_${count}_dir${dir}.png -c gray -m -15.5 -M -10.
map2img.py dens_${count}_dir${dir}.map -o dens_${count}_dir${dir}.png -l -c gray -m 0.15 -M 15000.
map2img.py temp_${count}_dir${dir}.map -o temp_${count}_dir${dir}.png -l -c gray -m 1e-5 -M 5e-3
map2img.py metal_${count}_dir${dir}.map -o metal_${count}_dir${dir}.png -l -c gray -m 2e-5 -M 2e-2

