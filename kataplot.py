#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import xarray
import pylab
import sys


datfile = sys.argv[1]

stokes = 'IQUV'
stokes = stokes.upper()

datfiles = []
for param in stokes:
	datfiles.append(datfile.replace('_I','_'+param))

for datfile in datfiles:

	data = xarray.open_dataset(datfile)
	corr = list(data.keys())[0]
	print('Plotting Stokes %s' % corr)

	pngname = datfile+'.png'

	fig = pylab.figure(figsize=(32,14))
	data[corr].plot.imshow(cmap='gist_heat',
							vmin = 0.0,
							vmax = 2e-2,
#							ylim = (0.55e9,1.02e9),
							interpolation='bicubic')

	fig.savefig(pngname,bbox_inches='tight')
