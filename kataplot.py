#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy
import os
import sys
import xarray
from optparse import OptionParser


def get_files(datfile):
	file_patterns = ['StokesI','StokesQ','StokesU','StokesV']
	to_plot = []
	for fp in file_patterns:
		if fp in datfile:
			input_type = fp
			file_patterns.remove(fp)
			to_plot.append(datfile)
	for fp in file_patterns:
		temp = datfile.replace(input_type,fp)
		if os.path.isfile(temp):
			to_plot.append(temp)
	if to_plot == []:
		to_plot.append(datfile)
		print('Cannot find other Stokes parameters, perhaps filename is non-standard?')
	return to_plot


def set_fontsize(fig,fontsize):
	def match(artist):
		return artist.__module__ == 'matplotlib.text'
	for textobj in fig.findobj(match=match):
		textobj.set_fontsize(fontsize)


# Set up command line options
parser = OptionParser(usage = '%prog [options] filename')
parser.add_option('--nostokes', dest = 'nostokes', help = 'Do not attempt to find and plot all Stokes, only the supplied file', action = 'store_true', default = False)
parser.add_option('--cmap', dest = 'mycmap', help = 'Colour map to use (default = gist_heat)', default = 'gist_heat')
parser.add_option('--flagcol', dest = 'flagcol', help = 'Colour to use for flagged data (default = black)', default = 'black')
parser.add_option('--interp', dest = 'interp', help = 'Image interpolation (default = bicubic)', default = 'bicubic')
parser.add_option('--vmin', dest = 'vmin', help = 'Minimum value for colour scale (default = 0.0)', default = 0.0)
parser.add_option('--vmax', dest = 'vmax', help = 'Maximum value for colour scale (default = half of peak value)', default = '')
parser.add_option('--tmin', dest = 'tmin', help = 'Minimum of time range (data units, default = minimum)', default = '')
parser.add_option('--tmax', dest = 'tmax', help = 'Maximum of time range (data units, default = maximum)', default = '')
parser.add_option('--fmin', dest = 'fmin', help = 'Minimum of frequency range (data units, default = minimum)', default = '')
parser.add_option('--fmax', dest = 'fmax', help = 'Maximum of frequency range (data units, default = maximum)', default = '')
parser.add_option('--xsize', dest = 'xsize', help = 'Figure width (default = 32)', default = 32)
parser.add_option('--ysize', dest = 'ysize', help = 'Figure width (default = 16)', default = 16)
parser.add_option('--fontsize', dest = 'fontsize', help = 'Font size for all figure elements (default = 20)', default = 20)


# Get command line options
(options,args) = parser.parse_args()
nostokes = options.nostokes
mycmap = options.mycmap
flagcol = options.flagcol
interp = options.interp
vmin = float(options.vmin)
vmax = options.vmax
tmin = options.tmin
tmax = options.tmax
fmin = options.fmin
fmax = options.fmax
xsize = float(options.xsize)
ysize = float(options.ysize)
fontsize = options.fontsize


# Check that the essential one is there
if len(args) != 1:
    print('Please specify a file to plot')
    sys.exit()
else:
    datfile = args[0]


# Find all the other Stokes params (or not)
if not nostokes:
	to_plot = get_files(datfile)


# Set up colour map and NaN colour
cmap = mpl.cm.get_cmap(mycmap).copy()
cmap.set_bad(flagcol,1.0)


# Loop over files
for datfile in to_plot:

	data = xarray.open_dataset(datfile)
	corr = list(data.keys())[0]
	pngname = datfile+'.png'
	print('Plotting Stokes %s (%s)' % (corr,datfile)) 

	# Sort out vmax if not supplied by user
	if vmax != '':
		vmax = float(vmax)
	else:
		raw = data[corr].values
		vmax = 0.5*numpy.max(raw[~numpy.isnan(raw)])

	# Frequency plot limits
	if fmin == '':
		fmin = numpy.min(data[corr].freq.values)
	else:
		fmin = float(fmin)
	if fmax == '':
		fmax = numpy.max(data[corr].freq.values)
	else:
		fmax = float(fmax)

	# Time plot limits
	if tmin == '':
		tmin = numpy.min(data[corr].time.values)
	else:
		tmin = float(tmin)
	if tmax == '':
		tmax = numpy.max(data[corr].time.values)
	else:
		tmax = float(tmax)


	fig = plt.figure(figsize=(xsize,ysize))
	data[corr].plot.imshow(cmap=cmap,
							vmin = vmin,
							vmax = vmax,
							xlim = (tmin,tmax),
							ylim = (fmin,fmax),
							interpolation=interp)

	set_fontsize(fig,fontsize)
	fig.savefig(pngname,bbox_inches='tight')
