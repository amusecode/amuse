#!/bin/sh

if (test -z $1); then
   echo Syntax is:
   echo    $0 name1 [name2 [name3 [...] ] ]
   exit 1;
fi

lastname=""
echo "Doing stages "$*
for name in $*; do
   echo "   Evolving stage "$name
   cp $name.dat init.dat
   cp $name.run init.run
   
   if (test -e $lastname.last1); then
      cp $lastname.last1 $name.start
   fi
   
   lastname=$name
   
   ./continue.sh $name
done
echo Stages done.
