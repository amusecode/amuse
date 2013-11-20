#!/usr/bin/python
# Draw multigrid levels as dumped by dump_mg_levels
 
import sys
from pyx import *

infile = sys.argv[1]

f = file(infile,'r')

scl = 20.0

ndim   = int(f.readline())
myid   = int(f.readline())
ncpu   = int(f.readline())
ilevel = int(f.readline())

if (ndim != 2):
	raise "Only works for ndim=2!"

print 'ndim   = ',ndim
print 'myid   = ',myid
print 'myid   = ',myid,'/',ncpu
print 'ilevel = ',ilevel

# Loop over levels
for i in range(1,ilevel):
	c = canvas.canvas()
	c.stroke(path.rect(0,0,scl,scl),[color.rgb.red])

	dx = 0.5**(i-1)*scl
	# Read active grids
	ngrids = int(f.readline())
	for igrid in range(ngrids):
		x = float(f.readline())*scl
		y = float(f.readline())*scl
		pt = path.rect(x-dx*0.5, y-dx*0.5, dx, dx)
		c.stroke(pt,[color.rgb.black, deco.filled([color.rgb(0.8,0.8,0.8)])])

	# Read reception grids
	for icpu in range(ncpu):
		if(icpu==myid-1): continue
		ngrids = int(f.readline())
		for igrid in range(ngrids):
			x = float(f.readline())*scl
			y = float(f.readline())*scl
			pt = path.rect(x-dx*0.5, y-dx*0.5, dx, dx)
			fcol = color.hsb(h=float(icpu)/ncpu,s=0.5,b=1.0)
			c.stroke(pt,[color.rgb.black, deco.filled([fcol])])
			c.text(x,y,str(icpu+1),[text.halign.boxcenter,text.valign.middle])

	c.writePDFfile(infile+'_level_'+str(i))

