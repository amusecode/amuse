#! /bin/tcsh
echo $WORKDIR
cd $WORKDIR
mkdir tmp_distrib
cd tmp_distrib

echo $ROOT
cp -r $ROOT ramses
cd ramses
rm -rf svn-commit*
rm -rf .svn
rm -rf */.svn
rm -rf */*/.svn
rm -rf */*/*/.svn
rm -rf */*/*/*/.svn

rm -rf *~
rm -rf */*~
rm -rf */*/*~
rm -rf */*/*/*~
rm -rf */*/*/*~/*~

rm -rf diffusion
rm -rf multimat
rm -rf stiff
rm -rf grafic2
rm -rf galic

cd patch
rm -rf diffusion
rm -rf induction

cd ../doc
cp src/ramses_ug.pdf .
rm -rf src

cd ../namelist
rm -rf *_mmat.nml
rm -rf *_stiff.nml
rm -rf *_diff.nml
rm -rf alloy.nml beltrami.nml ponomarenko.nml implosion.nml spitzer.nml

cd ../bin
rm -rf *.o
rm -rf *.mod
rm -rf *1d
rm -rf *2d
rm -rf *3d
rm -rf *_grafic
rm -rf histo
rm -rf map_gas*
rm -rf map_mass*
rm -rf amrdir
rm -rf sod
rm -rf amr2map
rm -rf amr2cube
rm -rf amr2gmsh
rm -rf part2map
rm -rf part2cube
rm -rf header
rm -rf log2col
mv Makefile.ramses Makefile
rm -rf Makefile.*

cd ../utils/f90
rm -rf *_grafic
rm -rf histo
rm -rf test_io
rm -rf *.o
rm -rf *.mod
rm -rf map_gas*
rm -rf map_mass*
rm -rf amr2cell
rm -rf amr2map
rm -rf amrdir
rm -rf sod
rm -rf amr2cube
rm -rf amr2gmsh
rm -rf part2map
rm -rf part2cube
rm -rf ramses2tipsy
rm -rf sunset
rm -rf vrot
rm -rf header
rm -rf log2col
rm -rf amr2cut
rm -rf amr2cylprof
rm -rf amr2prof
rm -rf amr2tipsy
rm -rf dbl2sng
rm -rf defrag
rm -rf getstarlist
rm -rf io_ramses
rm -rf part2cylprof
rm -rf part2prof
rm -rf part2ngpgrafic
rm -rf part2tipsy
rm -rf random
rm -rf hop_ramses/hop
rm -rf hop_ramses/regroup
rm -rf hop_ramses/poshalo
rm -rf hop_ramses/*.o
rm -rf zoom_ic/geticref
rm -rf zoom_ic/icdegrade
rm -rf zoom_ic/icinject
rm -rf zoom_ic/icupgrade

cd ../idl
rm -rf smooth

cd ../../..

tar cvf ramses.tar ramses
gzip ramses.tar

echo 'ramses package done' 
