#!/bin/sh
# Modified by Evert Glebbeek and converted to csh
# Allows to specify destination path for output files

# Get commandline arguments
if (test -z $2); then
   if (test -z $1); then
      echo Syntax is:
      echo    $0 [destpath] name
      echo where destpath is the destination directory \(default: current\) and name is the 
      echo output name
      exit 1;
   else
      outpath=`pwd`
      name=$1
   fi
else
      outpath=$1
      name=$2
fi

# Locate the evolution program
if (test -z $EV); then
   EV=`pwd`/bin/linux/ev
fi

RM=rm

# Remove files from previous run
if (test -e $name.out1); then 
   $RM $name.out1 
fi
if (test -e $name.out2); then
    $RM $name.out2
fi
if (test -e $name.log); then
    $RM $name.log
fi
if (test -e $name.out); then
    $RM $name.out
fi
if (test -e $name.last1); then
    $RM $name.last1
fi
if (test -e $name.last2); then
    $RM $name.last2
fi
if (test -e $name.mod); then
    $RM $name.mod
fi
if (test -e $name.plt1); then
    $RM $name.plt1
fi
if (test -e $name.plt2); then
    $RM $name.plt2
fi
if (test -e $name.mdl1); then
    $RM $name.mdl1
fi
if (test -e $name.mdl2); then
    $RM $name.mdl2
fi

# Create directories if needed
mkdir -p $outpath

# config files
ln -s init.dat fort.22
ln -s init.run fort.23

# input files
ln -s $outpath/input/zahb.mod fort.12
ln -s $outpath/input/zams.mod fort.16
ln -s $outpath/input/zams.dat fort.17
ln -s $outpath/input/zams.out fort.18
ln -s $outpath/input/zams.mas fort.19
ln -s $outpath/input/phys.z02 fort.20
ln -s $outpath/input/lt2ubv.dat fort.21

# output files
ln -s $name.out1 fort.1
ln -s $name.out2 fort.2
ln -s $name.io12 fort.3
ln -s $name.log fort.8
ln -s $name.out fort.9
ln -s $name.last1 fort.13
ln -s $name.last2 fort.14
ln -s $name.mod fort.15
ln -s $name.plt1 fort.31
ln -s $name.plt2 fort.32
ln -s $name.mdl1 fort.33
ln -s $name.mdl2 fort.34

# run code
#nice -n 19 $EV $name
$EV $name
#gdb $EV

# remove links
rm -f fort.16
rm -f fort.17
rm -f fort.18
rm -f fort.19
rm -f fort.20
rm -f fort.21
rm -f fort.22
rm -f fort.23

rm -f fort.1
rm -f fort.2
rm -f fort.3
rm -f fort.8
rm -f fort.9
rm -f fort.12
rm -f fort.13
rm -f fort.14
rm -f fort.15
rm -f fort.31
rm -f fort.32
rm -f fort.33
rm -f fort.34
