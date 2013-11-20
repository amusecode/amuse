#!/usr/bin/env python
 
import sys
import numpy
import pylab

import Image
import fortranfile

from optparse import OptionParser

def main():

	# Parse command line arguments
	parser = OptionParser()
	parser.usage = "%prog [options] map_file"
	parser.add_option('-l','--logscale',dest='logscale', action='store_true', \
	       help='use log color scaling', default=False)
	parser.add_option("-o","--output",  dest="outfile", metavar="FILE", \
			help='output image file [default: <map_file>.png]', default=None)
	parser.add_option("-m","--min",  dest="min", metavar="VALUE", \
			help='min value', default=None)
	parser.add_option("-M","--max",  dest="max", metavar="VALUE", \
			help='max value', default=None)
	parser.add_option('-a','--autorange',dest='autorange', action='store_true', \
	       help='use automatic dynamic range (overrides min & max)', default=False)
	parser.add_option('--big-endian',dest='big_endian', action='store_true', \
	       help='input binary data is stored as big endian', default=False)
	parser.add_option('-c','--colormap',dest='cmap_str', metavar='CMAP', \
	       help='matplotlib color map to use', default="jet")
	(opts,args)=parser.parse_args()

	# Parse input and output
	try:
		infile  = args[0]
	except:
		print parser.print_help()
		return 1

	if(opts.outfile==None):
		outfile=infile+'.tif'
	else:
		outfile=opts.outfile

	# Endianness
	if(opts.big_endian):
		endianness = ">"
	else:
		endianness = "="

	# Read image data
	print "Reading raw Fortran data..."
	f = fortranfile.FortranFile(infile)
	[t,dx,dy,dz] = f.read_fortran_record('f8', endian=endianness)
	[nx,ny] = f.read_fortran_record('i4', endian=endianness)
	dat = f.read_fortran_record('f4', endian=endianness)
	f.close()

	if(opts.logscale):
		dat = numpy.array(dat)+1e-12

	rawmin = numpy.amin(dat)
	rawmax = numpy.amax(dat)
	print '    Image map size  : ',(nx, ny)
	print '    Data bounds     : ',(rawmin,rawmax)

	print "Scaling data and processing colormap..."

	# Bounds
	if opts.min==None:
		plotmin = rawmin
	else:
		plotmin = float(opts.min)

	if opts.max==None:
		plotmax = rawmax
	else:
		plotmax = float(opts.max)

	# Log scale?
	if(opts.logscale):
		dat = numpy.log10(dat)
		rawmin = numpy.log10(rawmin)
		rawmax = numpy.log10(rawmax)
		plotmin = numpy.log10(plotmin)
		plotmax = numpy.log10(plotmax)

	# Auto-adjust dynamic range?
	if(opts.autorange):
		print "Computing dynamic range..."
		# Overrides any provided bounds
		NBINS = 200
		# Compute histogram
		(hist,bins) = numpy.histogram(dat, NBINS, (rawmin,rawmax), normed=True)
		chist = numpy.cumsum(hist); chist = chist / numpy.amax(chist)
		# Compute black and white point
		clip_k = chist.searchsorted(0.05)
		plotmin = bins[clip_k]
		plotmax = rawmax

	if(plotmax-plotmin>0):
		dat = numpy.clip((dat-plotmin)/(plotmax-plotmin), 0.0, 1.0)
	else:
		dat = 0.5*dat/plotmax

	if(opts.logscale):
		print '    Color bounds    : ',(10**plotmin,10**plotmax)
	else:
		print '    Color bounds    : ',(plotmin,plotmax)

	# Apply chosen color map
	color_map = pylab.cm.get_cmap(opts.cmap_str)
	dat = color_map(dat)*255

	# Convert to int
	dat = numpy.array(dat, dtype='i')

	# Output to file
	print "Saving image to file..."
	R_band = Image.new("L",(nx,ny))
	R_band.putdata(dat[:,0])
	G_band = Image.new("L",(nx,ny))
	G_band.putdata(dat[:,1])
	B_band = Image.new("L",(nx,ny))
	B_band.putdata(dat[:,2])

	out_img = Image.merge("RGB", (R_band, G_band, B_band)).transpose(Image.FLIP_TOP_BOTTOM)
	out_img.save(outfile)

if __name__ == '__main__':
	main()
