#Build galaxy models
BIN=../../bin


galaxy: disk1 halo1 bulge1
	cat disk1 bulge1 halo1 > galaxy

checkvirial: galaxy
	$(BIN)/checkvirial dbh.dat < galaxy

potential: dbh.dat

galaxy_ascii: galaxy
	$(BIN)/toascii galaxy > galaxy_ascii.dat

disk: freqdbh.dat cordbh.dat dbh.dat in.disk
	$(BIN)/gendisk < in.disk > disk

bulge: dbh.dat in.bulge
	$(BIN)/genbulge < in.bulge > bulge 

halo: dbh.dat in.halo
	$(BIN)/genhalo < in.halo > halo

disk1: disk dbh.dat
	$(BIN)/energysort dbh.dat disk disk > disk1

bulge1: bulge dbh.dat
	$(BIN)/energysort dbh.dat bulge bulge > bulge1

halo1: halo dbh.dat
	$(BIN)/energysort dbh.dat halo halo > halo1

dbh.dat: in.dbh
	$(BIN)/dbh < in.dbh

freqdbh.dat: dbh.dat h.dat b.dat
	$(BIN)/getfreqs

cordbh.dat: dbh.dat freqdbh.dat in.diskdf
	$(BIN)/diskdf < in.diskdf

vr.dat: dbh.dat h.dat b.dat
	$(BIN)/vcirc 1 1 > vr.dat

clean:
	rm -f disk bulge halo galaxy disk1 bulge1 halo1

veryclean:
	rm -f *.dat *.out disk bulge halo galaxy galaxy1 disk1 bulge1 halo1 toomre2.5 dfhalo.table
