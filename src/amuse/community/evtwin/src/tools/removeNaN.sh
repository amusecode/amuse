#! /bin/csh -f
###############################################################
#  removeNaN.sh
#  removes lines containing NaN (Not a Number) from .plt files
#  Selma de Mink 28 June 2005
#
###############################################################

# make outfile the same as infile to overwrite the infile
set infile = $1.plt$2
set outfile = $1.plt$2_withoutNaN


# check users input (never trust a human being)
# check script was called with two arguments
if ($#argv != 2) then
	echo "error: two areguments are expected"
	echo "call with: ./removeNaN.sh arg1 arg2"
	echo "    arg1 = <plt-filename without .plt>"
	echo "    arg2 = <1 or 2>, 1 for plt1, etc"
else

# check if the file is really a file
if (! -f $infile) then
	echo "error: $infile is not a file" 
else	
# check if the file exists
if (! -e $infile) then
	echo "error: $infile does not exist" 
else

# check silently if file contains NaN
echo "... checking file" $1.plt$2 " ..."
grep -s "NaN" $1.plt$2 >> /dev/null

if ($status == 0) then
	echo "... NaN was found in your file ..."
	grep -v "NaN" $infile > $outfile
	echo "... "$outfile " contains clean file ..."
else
	echo "No NaNs found in your file"
endif
endif
endif
endif
