#for ((a=10; a <= 99 ; a=a+1))
#do
#  echo output_000$a
#  hop -in output_000$a/part_000$a.out -p 1. -o hop000$a
#  regroup -root hop000$a -douter 80. -dsaddle 200. -dpeak 240. -f77 -o grp000$a
#  poshalo -inp output_000$a -pre grp000$a -xmi 0.375 -xma 0.625 -ymi 0.375 -yma 0.625 -zmi 0.375 -zma 0.625
#done

for ((a=140; a <= 347 ; a=a+1))
do
  echo output_00$a
  hop -in output_00$a/part_00$a.out -p 1. -o hop00$a
  regroup -root hop00$a -douter 80. -dsaddle 200. -dpeak 240. -f77 -o grp00$a
  poshalo -inp output_00$a -pre grp00$a -xmi 0.375 -xma 0.625 -ymi 0.375 -yma 0.625 -zmi 0.375 -zma 0.625
done


