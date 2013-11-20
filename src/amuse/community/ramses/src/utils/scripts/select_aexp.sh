#!/bin/bash
find $1 -name 'info*txt' -exec grep -H aexp '{}' + > dummy
sort dummy -o dummy_sorted
rm dummy
step=`echo $2 | awk '{print 1./$1}'`
asel_old=0.0
for ((n=1; n<=$2; n+=1)); do
    asel=`echo $asel_old $step | awk '{print $1+$2}'`
    mindiff=1000.
    mindir='xxx'
    while read line; do 
       tmp=`echo $line | cut -d = -f 1`
       len=`echo ${#tmp}`
# remove 21 chars at end because that is exactly the string '/info_NNNNN.txt:aexp '
       rest=`echo $len | awk '{print $1-21}'`
       curdir=`expr substr $tmp 1 $rest`
       anow=`echo $line | cut -d = -f 2`
       adiff=`echo $anow $asel | awk '{print $1-$2}'| awk ' { if($1>=0) {print $1} else {print $1*-1 }}'`
       s=`   echo $mindiff $adiff                 | awk '{if ($1 < $2) {print $1} else {print $2}}'`
       sdir=`echo $mindiff $adiff $mindir $curdir | awk '{if ($1 < $2) {print $3} else {print $4}}'`
       mindiff=$s
       mindir=$sdir
    done < dummy_sorted
 
    asel_old=$asel
    
    echo $mindir
done
rm dummy_sorted