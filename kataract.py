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
    logging.info('MS has %s channels' % len(chan_freqs))
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
    logging.info('MS has %s timeslots' % len(times))
    return times


def regularise_times(times):
    new_times = []
    time_flags = []
    diffs = [times[i + 1] - times[i] for i in range(len(times) - 1)]
    t_int = numpy.min(numpy.unique(numpy.round(diffs, 3)))
    for i in range(0,len(times)-1):
        new_times.append(times[i])
        time_flags.append(False)
        if numpy.round((times[i+1] - times[i]),3) != t_int:
            num_gaps = int(((times[i+1] - times[i]) // t_int))
            for j in range(1,num_gaps+1):
#                times.insert(i+j,times[i]+(j*t_int))
                new_times.append(times[i]+(j*t_int))
                time_flags.append(True)
    new_times.append(times[-1])
    time_flags.append(False)
    return new_times,time_flags


def save_xarray(datacube,name,chan_freqs,times,opfile):
    '''
    Put per-Stokes data into xarray and dump to netCDF file
    Rotate input array so time moves left to right when plotted
    '''
    logging.info('Writing %s' % opfile)
    datacube = numpy.stack(datacube)
    datacube = numpy.swapaxes(datacube,0,1)
    xdatacube = xarray.DataArray(datacube,
        name = name,
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
    parser.add_option('--minuv', dest = 'minuv', help = 'Minimum UV distance (in metres) to use (default = 0)', default = 0)
    parser.add_option('--t0', dest = 't0', help = 'First time slot to plot (default = plot all)', default = 0)
    parser.add_option('--t1', dest = 't1', help = 'Final time slot to plot (default = plot all)', default = 0)
    parser.add_option('--subtract-col', dest = 'subtract_col', help = 'Subtract this column on the fly (default = no subtraction)', default = False)
    parser.add_option('--ignore-weights', dest = 'ignore_weights', help = 'Ignore WEIGHT_SPECTRUM when computing spectra (default = False)', action = 'store_true', default = False)
    parser.add_option('--flag-threshold', dest = 'flag_threshold', help = 'Percentage of flags per time/freq point above which datum is discarded (default = 50)', default = 50)
    parser.add_option('--ifile', dest = 'ifile', help = 'netCDF file for Stokes I output (default = based on MS name)', default = '')
    parser.add_option('--qfile', dest = 'qfile', help = 'netCDF file for Stokes Q output (default = based on MS name)', default = '')
    parser.add_option('--ufile', dest = 'ufile', help = 'netCDF file for Stokes U output (default = based on MS name)', default = '')
    parser.add_option('--vfile', dest = 'vfile', help = 'netCDF file for Stokes V output (default = based on MS name)', default = '')
    parser.add_option('--xxfile', dest = 'xxfile', help = 'netCDF file for XX output (default = based on MS name)', default = '')
    parser.add_option('--xyfile', dest = 'xyfile', help = 'netCDF file for XY output (default = based on MS name)', default = '')
    parser.add_option('--yxfile', dest = 'yxfile', help = 'netCDF file for YX output (default = based on MS name)', default = '')
    parser.add_option('--yyfile', dest = 'yyfile', help = 'netCDF file for YY output (default = based on MS name)', default = '')
    parser.add_option('--nocorr', dest = 'nocorr', help = 'Do not save baseline-averaged correlations (default = save them)', default = False, action = 'store_true')
    parser.add_option('--logfile', dest = 'logfile', help = 'File name for log output (default = based on MS name)', default = '')
    parser.add_option('--init-log', dest = 'initlogfile', help = 'Remove existing log file at the start (default = do not remove)', action = 'store_true', default = False)

    # Get command line options
    (options,args) = parser.parse_args()
    mycol = options.mycol
    stokes = options.stokes.upper()
    myfield = int(options.myfield)
    myspw = int(options.myspw)
    minuv = float(options.minuv)
    t0 = int(options.t0)
    t1 = int(options.t1)
    subtract_col = options.subtract_col
    ignore_weights = options.ignore_weights
    flag_threshold = options.flag_threshold
    ifile = options.ifile
    qfile = options.qfile
    ufile = options.ufile
    vfile = options.vfile
    xxfile = options.xxfile
    xyfile = options.xyfile
    yxfile = options.yxfile
    yyfile = options.yyfile
    nocorr = options.nocorr
    logfile = options.logfile
    initlogfile = options.initlogfile

    # Check that the essential one is there
    if len(args) != 1:
        print('Please specify a Measurement Set')
        sys.exit()
    else:
        myms = args[0].rstrip('/')

    # Setup output file names
    field_label = '_field'+str(myfield)
    if ifile == '':
        ifile = myms.split('/')[-1]+field_label+'_StokesI-kata.nc'
    if qfile == '':
        qfile = myms.split('/')[-1]+field_label+'_StokesQ-kata.nc'
    if ufile == '':
        ufile = myms.split('/')[-1]+field_label+'_StokesU-kata.nc'
    if vfile == '':
        vfile = myms.split('/')[-1]+field_label+'_StokesV-kata.nc'
    if xxfile == '':
        xxrfile = myms.split('/')[-1]+field_label+'_XX-kata.nc'
    if xyfile == '':
        xyrfile = myms.split('/')[-1]+field_label+'_XY-kata.nc'
    if yxfile == '':
        yxrfile = myms.split('/')[-1]+field_label+'_YX-kata.nc'
    if yyfile == '':
        yyrfile = myms.split('/')[-1]+field_label+'_YY-kata.nc'

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
    logging.info('      ====== kataract ======')
    logging.info('')
    logging.info('Processing %s ' % myms)

    # Make list of required columns for daskms
    logging.info('Visibilities will be read from the %s column' % mycol)
    ms_cols = [mycol,'FLAG']
    if subtract_col:
        logging.info('Residuals formed by subtraction of %s column' % subtract_col)
        ms_cols.append(subtract_col)
    if minuv > 0.0:
        logging.info('Baselines below %s m will be excluded' % minuv)
    if not ignore_weights:
        ms_cols.append('WEIGHT_SPECTRUM')
    else:
        logging.info('WEIGHT_SPECTRUM will be ignored')

    logging.info('Getting MS info, please wait...')
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

    if not nocorr:
        xxcube =[]
        xycube =[]
        yxcube =[]
        yycube =[]

    logging.info('Getting data...')
    # Loop over timeslots
    start_time = time.time()
    for i in range(len(times)):

        timeslot = times[i]
        # Get data and flags for this timeslot
        msdata = daskms.xds_from_ms(myms, columns = ms_cols, 
            taql_where = 'TIME_CENTROID==%s && DATA_DESC_ID==%s  && FIELD_ID==%s && sqrt(sumsqr(UVW[:2])) > %s' % (timeslot,myspw,myfield,minuv))
        
        vis = msdata[myspw][mycol].values
        flags = msdata[myspw]['FLAG'].values
        if not ignore_weights:
            weights = msdata[myspw]['WEIGHT_SPECTRUM'].values

        # Parallel reads for a zarr-backed MS
        # vis = msdata[myspw][mycol].data
        # flags = msdata[myspw]['FLAG'].data
        # weights = msdata[myspw]['WEIGHT_SPECTRUM'].data
        # vis,flags,weights = dask.compute(vis,flags,weights)

        # Collapse flags across baseline axis and turn into a percentage
        # Resulting array has shape (nchan,ncorr)
        flag_mask = (100.0*numpy.mean(flags,axis=0) > flag_threshold)

        # Compute residuals if requested
        if subtract_col:
            vis = vis - msdata.MODEL_DATA.values
        
        # Apply flags as a mask and compute per-corr weighted vector averaged amplitudes
        vis = numpy.ma.array(vis, mask = flags)
        # XX = numpy.ma.mean(vis[:,:,0],axis=0)
        # XY = numpy.ma.mean(vis[:,:,1],axis=0)
        # YX = numpy.ma.mean(vis[:,:,2],axis=0)
        # YY = numpy.ma.mean(vis[:,:,3],axis=0)
        if not ignore_weights:
            weights = numpy.ma.array(weights, mask = flags)
            XX = numpy.ma.average(vis[:,:,0], axis=0, weights=weights[:,:,0])
            XY = numpy.ma.average(vis[:,:,1], axis=0, weights=weights[:,:,1])
            YX = numpy.ma.average(vis[:,:,2], axis=0, weights=weights[:,:,2])
            YY = numpy.ma.average(vis[:,:,3], axis=0, weights=weights[:,:,3])
        else:
            XX = numpy.ma.average(vis[:,:,0], axis=0)
            XY = numpy.ma.average(vis[:,:,1], axis=0)
            YX = numpy.ma.average(vis[:,:,2], axis=0)
            YY = numpy.ma.average(vis[:,:,3], axis=0)

        # Apply flag percentage cuts
        XX.data[numpy.where(flag_mask[:,0]==True)] = numpy.nan
        XY.data[numpy.where(flag_mask[:,1]==True)] = numpy.nan
        YX.data[numpy.where(flag_mask[:,2]==True)] = numpy.nan
        YY.data[numpy.where(flag_mask[:,3]==True)] = numpy.nan

        xxcube.append(XX)
        xycube.append(XY)
        yxcube.append(YX)
        yycube.append(YY)

        ispectrum = numpy.real((XX+YY)/2.0)
        qspectrum = numpy.real((XX-YY)/2.0)
        uspectrum = numpy.real((XY+YX)/2.0)
        vspectrum = numpy.imag((XY-YX)/(2.0*j))

        if 'I' in stokes:
            icube.append(ispectrum)
        if 'Q' in stokes:
            qcube.append(qspectrum)
        if 'U' in stokes:
            ucube.append(uspectrum)
        if 'V' in stokes:
            vcube.append(vspectrum)

        count += 1
        if count == ten_percent:
            end_time = time.time()
            elapsed = round(end_time-start_time,1)
            eta = round(9*elapsed,1)
            logging.info('First ten percent of spectra took %s seconds' % elapsed)
            logging.info('Estimated completion in %s seconds' % eta)

    
    logging.info('Regularising time axis')
    all_times,time_flags = regularise_times(times)
    dummy_spectrum = numpy.ones(len(chan_freqs))*numpy.nan

    icube_final = []
    qcube_final = []
    ucube_final = []
    vcube_final = []

    j = 0
    for i in range(0,len(time_flags)):
        if time_flags[i]:
            icube_final.append(dummy_spectrum)
            qcube_final.append(dummy_spectrum)
            ucube_final.append(dummy_spectrum)
            vcube_final.append(dummy_spectrum)
        else:
            icube_final.append(icube[j])
            qcube_final.append(qcube[j])
            ucube_final.append(ucube[j])
            vcube_final.append(vcube[j])
            j+=1 

    times = all_times

    # Dump dynamic spectrum data to xarrays and save them
    logging.info('Saving requested xarrays')
    if 'I' in stokes:
        save_xarray(icube_final,'I',chan_freqs,times,ifile)
    if 'Q' in stokes:
        save_xarray(qcube_final,'Q',chan_freqs,times,qfile)
    if 'U' in stokes:
        save_xarray(ucube_final,'U',chan_freqs,times,ufile)
    if 'V' in stokes:
        save_xarray(vcube_final,'V',chan_freqs,times,vfile)
    if not nocorr:
        logging.info('Correlation dumps not implemented yet (netCDF has no complex support)')
        # save_xarray(xxcube,'XX',chan_freqs,times,xxfile)
        # save_xarray(xycube,'XY',chan_freqs,times,xyfile)
        # save_xarray(yxcube,'YX',chan_freqs,times,yxfile)
        # save_xarray(yycube,'YY',chan_freqs,times,xxfile)


    logging.info('Finished')

