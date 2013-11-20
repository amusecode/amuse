a=$1
if [ $a -ge 10 ]; then
    count='000'$a
fi
if [ $a -ge 100 ]; then
    count='00'$a
fi
if [ $a -ge 1000 ]; then
    count='0'$a
fi
echo output_$count

echo "'./output_$count/'   Ra3     1      $a" > inputfiles_HaloMaker.dat

source params_$count.sh

if [ $2 -eq 1 ]; then
    ~/HaloMaker_stars/HaloMaker ./
    mv tree_bricks tree_bricks_$count
fi

~/HaloMaker_stars/read_halo -inp output_$count -xmi $xc -ymi $yc -zmi $zc -rma $rvir > sat_$count.list
sed -i '1,8d' sat_$count.list

cat sat_$count.list







