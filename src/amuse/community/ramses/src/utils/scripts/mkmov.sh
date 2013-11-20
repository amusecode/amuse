#! /bin/csh -f
   set target = $argv[1]
   echo $target
   mkdir $target
   set wordlist = output_*
   foreach i ( $wordlist ) 
	echo $i
	amr2map -inp $i -out $target/map_$i.dat -typ 1
   end
