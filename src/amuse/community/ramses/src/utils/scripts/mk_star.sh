for ((n=1; n<124; n++ ))
  do 
  char=`sed -n ${n}p time_list.txt | cut -c 1-12`
  echo ${n}' '${char}
  $HOME/dom/ramses/utils/f90/amr2cube -inp ${char} -out cubes/dens_${char}.grafic -typ 1 -fil grafic -lma 10
done
