#!/usr/bin/env python
# ian.heywood@physics.ox.ac.uk


import daskms
import logging
import numpy
import os
import sys
import time
import xarray
from optparse import OptionParser


j = complex(0,1)



def get_chan_freqs(myms):
    '''
    Get a list containing channel frequencies
    '''
    spw_tab = daskms.xds_from_table(myms+'::SPECTRAL_WINDOW', columns = 'CHAN_FREQ')
    chan_freqs = spw_tab[0].CHAN_FREQ[myspw,:].values
    logging.info('%s channels' % len(chan_freqs))
    return chan_freqs


def get_times(myms,t0,t1):
    '''
    Get a list of the unique values in the TIME_CENTROID column of the input MS
    Trim to indices between t0 and t1 if necessary
    '''
    time_centroid = daskms.xds_from_ms(myms, columns=['TIME_CENTROID'],taql_where = 'FIELD_ID==%s && DATA_DESC_ID==%s' % (myfield,myspw))[0]
    times, indices = numpy.unique(time_centroid.TIME_CENTROID.values, return_index = True)
    # time_centroid = time_centroid.isel(row=indices)
    # times = time_centroid.TIME_CENTROID.values
    if t1 == 0: t1 = len(times)
    times = times[t0:t1]
    logging.info('%s unique timeslots' % len(times))
    return times


def save_xarray(datacube,chan_freqs,times,opfile):
    '''
    Put per-Stokes data into xarray and dump to netCDF file
    Rotate input array so time moves left to right when plotted
    '''
    datacube = numpy.stack(datacube)
    datacube = numpy.swapaxes(datacube,0,1)
    xdatacube = xarray.DataArray(datacube,
        name = 'I',
        dims = ['freq','time'],
        coords = {'freq' : chan_freqs,
                  'time' : times})
    xdatacube.to_netcdf(opfile)


if __name__ == '__main__':
    
    parser = OptionParser(usage = '%prog [options] ms')
    parser.add_option('--column', dest = 'mycol', help = 'Data column to plot (default = CORRECTED_DATA)', default = 'CORRECTED_DATA')
    parser.add_option('--stokes', dest = 'stokes', help = 'Stokes parameters to plot (default = IQUV)', default = 'IQUV')
    parser.add_option('--field', dest = 'myfield', help = 'FIELD_ID selection (default = 0)', default = 0)
    parser.add_option('--spw', dest = 'myspw', help = 'SPW ID selection (default = 0)', default = 0)
    parser.add_option('--t0', dest = 't0', help = 'First time slot to plot (default = plot all)', default = 0)
    parser.add_option('--t1', dest = 't1', help = 'Final time slot to plot (default = plot all)', default = 0)
    parser.add_option('--subtract-col', dest = 'subtract_col', help = 'Subtract this column on the fly (default = no subtraction)', default = False)
    parser.add_option('--ignore-weights', dest = 'ignore_weights', help = 'Ignore WEIGHT_SPECTRUM when computing spectra (default = False)', action = 'store_true', default = False)
    parser.add_option('--ifile', dest = 'ifile', help = 'netCDF file for Stokes I output (default = based on MS name)', default = '')
    parser.add_option('--qfile', dest = 'qfile', help = 'netCDF file for Stokes Q output (default = based on MS name)', default = '')
    parser.add_option('--ufile', dest = 'ufile', help = 'netCDF file for Stokes U output (default = based on MS name)', default = '')
    parser.add_option('--vfile', dest = 'vfile', help = 'netCDF file for Stokes V output (default = based on MS name)', default = '')
    parser.add_option('--logfile', dest = 'logfile', help = 'File name for log output (default = based on MS name)', default = '')
    parser.add_option('--init-log', dest = 'initlogfile', help = 'Remove existing log file at the start (default = do not remove)', default = False)

    # Get command line options
    (options,args) = parser.parse_args()
    mycol = options.mycol
    stokes = options.stokes.upper()
    myfield = int(options.myfield)
    myspw = int(options.myspw)
    t0 = int(options.t0)
    t1 = int(options.t1)
    subtract_col = options.subtract_col
    ignore_weights = options.ignore_weights
    ifile = options.ifile
    qfile = options.qfile
    ufile = options.ufile
    vfile = options.vfile
    logfile = options.logfile
    initlogfile = options.initlogfile

    # Check that the essential one is there
    if len(args) != 1:
        print('Please specify a Measurement Set')
        sys.exit()
    else:
        myms = args[0].rstrip('/')

    # Setup output file names
    if ifile == '':
        ifile = 'kataract_'+myms.split('/')[-1]+'_I.nc'
    if qfile == '':
        ifile = 'kataract_'+myms.split('/')[-1]+'_Q.nc'
    if ufile == '':
        ifile = 'kataract_'+myms.split('/')[-1]+'_U.nc'
    if vfile == '':
        ifile = 'kataract_'+myms.split('/')[-1]+'_V.nc'

    # Setup logfile name
    if logfile == '':
        logfile = myms.split('/')[-1]+'-kataract.log'

    # Remove existing logfile if needed
    if initlogfile and os.path.isfile(logfile):
        os.remove(logfile)

    # Set up logger
    logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s |  %(message)s', datefmt='%d/%m/%Y %H:%M:%S ')
    logging.getLogger().addHandler(logging.StreamHandler())

    logging.info('')
    logging.info('k a t a r a c t')
    logging.info('')
    logging.info('TO DO: Implement weights')
    logging.info('TO DO: Implement flag propagation')
    logging.info('')
    logging.info('Processing %s ' % myms)

    # Make like of required columns for daskms
    ms_cols = [mycol,'FLAG']
    if subtract_col:
        ms_cols.append(subtract_col)
    if not ignore_weights:
        ms_cols.append('WEIGHT_SPECTRUM')

    # Get channel frequencies
    chan_freqs = get_chan_freqs(myms)

    # Get timeslots, work out 10% for the time estimate
    times = get_times(myms,t0,t1)
    ten_percent = int(len(times)/10.0)
    count = 0

    # Initialise empty arrays to receive the data
    icube = []
    qcube = []
    ucube = []
    vcube = []


    # Loop over timeslots
    start_time = time.time()
    for timeslot in times:

        # Get data and flags for this timeslot
        msdata = daskms.xds_from_ms(myms, columns = ms_cols, 
            taql_where = 'TIME_CENTROID==%s && DATA_DESC_ID==%s ' % (timeslot,myspw))
        vis = msdata[myspw][mycol].values
        flags = msdata[myspw]['FLAG'].values

        # Compute residuals if requested
        if subtract_col:
            vis = vis - msdata.MODEL_DATA.values
        
        # Apply flags as a mask and compute per-corr vector averaged amplitudes
        vis = numpy.ma.array(vis, mask = flags)
        XX = numpy.ma.mean(vis[:,:,0],axis=0)
        XY = numpy.ma.mean(vis[:,:,1],axis=0)
        YX = numpy.ma.mean(vis[:,:,2],axis=0)
        YY = numpy.ma.mean(vis[:,:,3],axis=0)

        # 
        if 'I' in stokes:
            ispectrum = numpy.abs((XX+YY)/2.0)
            icube.append(ispectrum)
        if 'Q' in stokes:
            qspectrum = numpy.abs((XX-YY)/2.0)
            qcube.append(qspectrum)
        if 'U' in stokes:
            uspectrum = numpy.abs((XY-YX)/2.0)
            ucube.append(uspectrum)
        if 'V' in stokes:
            vspectrum = numpy.abs((XY+YX)/(2.0*j))
            vcube.append(vspectrum)

        count += 1
        if count == ten_percent:
            end_time = time.time()
            elapsed = round(end_time-start_time,1)
            eta = round(9*elapsed,1)
            logging.info('First ten percent took %s seconds' % elapsed)
            logging.info('Estimated completion in %s seconds' % eta)


    # Dump dynamic spectrum data to xarrays and save them
    if 'I' in stokes:
        save_xarray(icube,chan_freqs,times,ifile)
    if 'Q' in stokes:
        save_xarray(qcube,chan_freqs,times,qfile)
    if 'U' in stokes:
        save_xarray(ucube,chan_freqs,times,ufile)
    if 'V' in stokes:
        save_xarray(vcube,chan_freqs,times,vfile)


    logging.info('Finished')

