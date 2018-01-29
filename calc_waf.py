#!/usr/bin/env python
"""
SVN INFO: $Id$
Filename:     calc_waf.py
Author:       Damien Irving, d.irving@student.unimelb.edu.au
Description:  Calculates the wave activity flux (waf)
Reference:    Takaya, K., Nakamura, H., 2001
              A formulation of a phase-independent wave-activity flux
	      for stationary and migratory quasi-geostrophic eddies
	      on zonally-varying basic flow. J. Atmos. Sci., 58, 608--627

Input:        u wind, u wind climatology, u wind var name,
              v wind, v wind climatology, v wind var name,
	      geopotential height, geopotential height climatology, height var name,
              out_wafx_%s_filename,   out_wafy_%s_filename,
              level, start_year, end_year, start_month, end_month, start_day,
              end_day, climatology_year

Output:       waf x component, waf y component (out_wafx_%s_filename -> %s will be replaced with individual year)


Updates | By | Description
--------+----+------------
12 October 2012 | Damien Irving | Initial version.
29 January 2018 | Arulalan.T    | Second version.

"""

__version__ = '$Revision$'  ## Change this for git (this is a svn relevant command)

### Import required modules ###

import optparse
from optparse import OptionParser
from datetime import datetime
import os, subprocess, sys
import cdms2, MV2
import numpy
import netCDF4
import regrid2
from regrid2 import Regridder


def time_axis_check(axis1, axis2):
    """Checks whether the time axes of the input files are the same"""

    start_time1 = axis1.asComponentTime()[0]
    start_time1 = str(start_time1)
    start_year1 = start_time1.split('-')[0]

    end_time1 = axis1.asComponentTime()[-1]
    end_time1 = str(end_time1)
    end_year1 = end_time1.split('-')[0]

    start_time2 = axis2.asComponentTime()[0]
    start_time2 = str(start_time2)
    start_year2 = start_time2.split('-')[0]

    end_time2 = axis2.asComponentTime()[-1]
    end_time2 = str(end_time2)
    end_year2 = end_time2.split('-')[0]

    if (start_year1 != start_year2 or len(axis1) != len(axis2)):
        sys.exit('Input files do not all have the same time axis')


def xy_axis_check(axis1, axis2):
    """Checks whether the lat or lon axes of the input files are the same"""

    if (len(axis1) != len(axis2)):
        sys.exit('Input files do not all have the same %s axis' % (axis1.id))


def read_climatology(fname, vname, nlat, nlon, lev, tdate):
    """Reads the climatology data, assuming a 12 step monthly file format"""

    infile = cdms2.open(fname)
    daily_climatology = infile(vname, time=tdate, level=lev, order='tyx', squeeze=1)
    infile.close()

    return daily_climatology


def write_binary(data, outfile):
    """Takes a numpy data array and writes a corresponding binary data file"""

    fileobj = open(outfile, mode='wb')
    # fwrite(fileobj, data.size, data)
    data.filled(0).flatten().tofile(fileobj)
    fileobj.close()


def read_binary(infile, dims):
    """Reads a binary data file and puts data in an array"""

    fileobj = open(infile, mode='rb')

    datatype = 'f'
    size = dims[0] * dims[1] * dims[2]
    # read_data = fread(fileobj, size, datatype)
    read_data = numpy.fromfile(fileobj, dtype=datatype)
    read_data = read_data.reshape(dims)
    fileobj.close()

    return read_data


def regrid(data, in_lat, in_lon, oldgrid):
    """Takes input data and regrids to a 2.5 degree global grid"""

    #tnf_xy_onelevel.run requires North to South latitude axis
    latitude = cdms2.createAxis(numpy.arange(90, -92.5, -2.5, 'f'), id='latitude')  #-90,92.5,2.5

    latitude.designateLatitude()
    latitude.units = 'degrees_north'
    latitude.long_name = 'Latitude'
    latitude.standard_name = 'latitude'
    latitude.axis = 'Y'

    longitude = cdms2.createAxis(numpy.arange(0, 360, 2.5, 'f'), id='longitude')
    longitude.designateLongitude()
    longitude.units = 'degrees_east'
    longitude.long_name = 'Longitude'
    longitude.standard_name = 'longitude'
    longitude.axis = 'X'
    longitude.modulo = 360.
    longitude.topology = 'circular'

    newgrid = cdms2.createRectGrid(latitude, longitude)
    regridfunc = Regridder(oldgrid, newgrid)
    new_data = regridfunc(data)

    return new_data, latitude, longitude


def write_netcdf(fname_out, waf_data_ns, time_axis, lat_axis, lon_axis, sourcefile_text, outvar):
    """Writes the output netcdf file"""

    outfile = netCDF4.Dataset(
        fname_out, 'w', format='NETCDF3_CLASSIC')  ## Found error using cdo on abyss with format='NETCDF4'

    # Global attributes #

    setattr(outfile, 'Title', 'Calculated wave activity flux from U wind, V wind and geopotential height')
    setattr(outfile, 'Contact', 'Damien Irving (d.irving@student.unimelb.edu.au)')
    setattr(outfile, 'Reference', 'Takaya & Nakamura (2001) J. Atmos. Sci., 58, 608-627')
    setattr(outfile, 'Sourcefiles', sourcefile_text)
    setattr(outfile, 'history', '%s: Created using %s, format=NETCDF3_CLASSIC' %
            (datetime.now().strftime("%a %b %d %H:%M:%S %Y"), sys.argv[0]))

    # Dimensions #

    outfile.createDimension('time', None)
    outfile.createDimension('latitude', len(lat_axis))
    outfile.createDimension('longitude', len(lon_axis))

    times = outfile.createVariable('time', 'f4', ('time',))
    lats = outfile.createVariable('latitude', 'f4', ('latitude',))
    lons = outfile.createVariable('longitude', 'f4', ('longitude',))

    for att_name in time_axis.attributes.keys():
        setattr(times, att_name, time_axis.attributes[att_name])
    for att_name in lat_axis.attributes.keys():
        setattr(lats, att_name, lat_axis.attributes[att_name])
    for att_name in lon_axis.attributes.keys():
        setattr(lons, att_name, lon_axis.attributes[att_name])

    # Variable #

    out_data = outfile.createVariable(outvar, 'f4', ('time', 'latitude', 'longitude',), fill_value=9.999e+20)
    setattr(out_data, 'standard_name', outvar)
    setattr(out_data, 'units', 'm2 s-2')
    setattr(out_data, 'long_name', 'Wave activity flux, %s component' % (outvar[-1]))
    setattr(out_data, 'missing_value', 9.999e+20)
    setattr(out_data, 'history', 'Calculated wave activity flux from U wind, V wind and geopotential height')

    waf_data_ns = waf_data_ns.astype(numpy.float32)

    #reshape output data so latitude axis is south to north
    times[:] = time_axis
    lats[:] = lat_axis[::-1]
    lons[:] = lon_axis
    out_data[:] = waf_data_ns[:, ::-1, :]

    outfile.close()


def apply_mask(input_data, input_clim_u, time_axis):
    """Apply mask, because WAF only valid where the climatological mean flow is westerly"""

    # Get the time axis info #

    time_data = time_axis.asComponentTime()
    months = []
    for ii in range(0, len(time_data)):
        months.append(int(str(time_data[ii]).split('-')[1]))

    # Expand input climatology to match size of input data (with correct month matching) #

    ntime, nlat, nlon = numpy.shape(input_data)
    input_clim_u_expanded = numpy.ma.ones([ntime, nlat, nlon]) * 9.999e+20
    for i in range(0, ntime):
        month_index = months[i]
        input_clim_u_expanded[i, :, :] = input_clim_u[month_index - 1, :, :]

    # Define the mask #

    mask_u = numpy.where(input_clim_u_expanded > 1.0, 0,
                         1)  # exclude points where climtological wind has easterly component
    mask_inf = numpy.isinf(input_data)  # exclude infinity values
    mask_nan = numpy.isnan(input_data)  # exclude NaN values

    mask = numpy.sum([mask_u, mask_inf, mask_nan], axis=0)

    # Apply the mask #

    input_data_masked = MV2.array(numpy.ma.masked_array(input_data, mask=mask))

    return input_data_masked


def read_data(fname, vname, lev, tdate):
    """Reads the data from a typical input file"""

    infile = cdms2.open(fname)
    data = infile(vname, time=tdate, level=lev, order='tyx', squeeze=1)
    time = data.getTime()
    lat = data.getLatitude()
    lon = data.getLongitude()
    infile.close()
    return data.squeeze(), time, lat, lon


def main(fname_u, fname_uclim, vname_u, fname_v, fname_vclim, vname_v, fname_zg, fname_zgclim, vname_zg, fname_wafx,
         fname_wafy, lev, smon, emon, sday, eday, year, cyear):
    """Run the program"""

    ### Read the input data ###
    dsdate = '%s-%s-%s' % (year, smon, sday)
    dedate = '%s-%s-%s' % (year, emon, eday)
    csdate = '%s-%s-%s' % (cyear, smon, sday)
    cedate = '%s-%s-%s' % (cyear, emon, eday)
    hPa = float(lev)
    ## u wind, v wind, zg ##

    data_u, time_u, lat_u, lon_u = read_data(fname_u, vname_u, lev=hPa, tdate=(dsdate, dedate))
    data_v, time_v, lat_v, lon_v = read_data(fname_v, vname_v, lev=hPa, tdate=(dsdate, dedate))
    data_zg, time_zg, lat_zg, lon_zg = read_data(fname_zg, vname_zg, lev=hPa, tdate=(dsdate, dedate))
    ### Check that the input data are all on the same coordinate axes ###

    ## Time ##

    time_axis_check(time_u, time_v)
    time_axis_check(time_u, time_zg)

    ## Latitude ##

    xy_axis_check(lat_u, lat_v)
    xy_axis_check(lat_u, lat_zg)

    ## Longitude ##

    xy_axis_check(lon_u, lon_v)
    xy_axis_check(lon_u, lon_zg)

    ## Find grid characteristics ##
    ntime, nlat, nlon = data_u.shape
    input_grid = data_u.getGrid()

    ## u wind, v wind, zg climatology ##

    data_uclim = read_climatology(fname_uclim, vname_u, nlat, nlon, lev=hPa, tdate=(csdate, cedate))
    data_vclim = read_climatology(fname_vclim, vname_v, nlat, nlon, lev=hPa, tdate=(csdate, cedate))
    data_zgclim = read_climatology(fname_zgclim, vname_zg, nlat, nlon, lev=hPa, tdate=(csdate, cedate))

    ### Regrid the data ###
    data_u.setAxisList([time_u, lat_u, lon_u])
    data_v.setAxisList([time_v, lat_v, lon_v])
    data_zg.setAxisList([time_zg, lat_zg, lon_zg])
    data_u_regrid, lat_regrid, lon_regrid = regrid(data_u, lat_u, lon_u, input_grid)
    data_uclim_regrid, lat_regrid, lon_regrid = regrid(data_uclim, lat_u, lon_u, input_grid)

    data_v_regrid, lat_regrid, lon_regrid = regrid(data_v, lat_u, lon_u, input_grid)
    data_vclim_regrid, lat_regrid, lon_regrid = regrid(data_vclim, lat_u, lon_u, input_grid)

    data_zg_regrid, lat_regrid, lon_regrid = regrid(data_zg, lat_u, lon_u, input_grid)
    data_zgclim_regrid, lat_regrid, lon_regrid = regrid(data_zgclim, lat_u, lon_u, input_grid)

    ### Create the binary data files for input into the Fortran wap script ###

    write_binary(data_u_regrid, 'u.bin')
    write_binary(data_uclim_regrid, 'uclim.bin')
    write_binary(data_v_regrid, 'v.bin')
    write_binary(data_vclim_regrid, 'vclim.bin')
    write_binary(data_zg_regrid, 'zg.bin')
    write_binary(data_zgclim_regrid, 'zgclim.bin')

    ### Calculate the wave activity flux ###

    fout = open('answers.txt', 'w')
    print >> fout, '%s' % hPa
    print >> fout, 'zg.bin'
    print >> fout, 'zgclim.bin'
    print >> fout, 'u.bin'
    print >> fout, 'uclim.bin'
    print >> fout, 'v.bin'
    print >> fout, 'vclim.bin'
    print >> fout, 'wafx.bin'
    print >> fout, 'wafy.bin'

    fout.close()

    subprocess.call("./tnf_xy_onelevel.run < answers.txt", shell=True)

    ### Write the output netCDF file ###

    wafx_data_ns = read_binary('wafx.bin', [ntime, 73, 144])
    wafy_data_ns = read_binary('wafy.bin', [ntime, 73, 144])

    wafx_data_ns_masked = apply_mask(wafx_data_ns, data_uclim_regrid, time_u)
    wafy_data_ns_masked = apply_mask(wafy_data_ns, data_uclim_regrid, time_u)

    sourcefile_text = '%s, %s, %s, %s, %s, %s' % (fname_u, fname_uclim, fname_v, fname_vclim, fname_zg, fname_zgclim)

    write_netcdf(fname_wafx, wafx_data_ns_masked, time_u, lat_regrid, lon_regrid, sourcefile_text, 'wafx')
    write_netcdf(fname_wafy, wafy_data_ns_masked, time_u, lat_regrid, lon_regrid, sourcefile_text, 'wafy')

    ### Clean up ###

    os.system("rm answers.txt zg.bin zgclim.bin u.bin uclim.bin v.bin vclim.bin wafx.bin wafy.bin")


if __name__ == '__main__':

    ### Help and manual information ###

    usage = """usage: %prog [options] {u_file} {u_clim_file} {u_name}
{v_file} {v_clim_file} {v_name} {zg_file} {zg_clim_file} {zg_name}
{output_wafx_file} {output_wafy_file} {hPa} {syear} {eyear} {smon} {emon} {sday} {edy} {cyear}"""
    parser = OptionParser(usage=usage)

    parser.add_option(
        "-M",
        "--manual",
        action="store_true",
        dest="manual",
        default=False,
        help="output a detailed description of the program")

    (options, args
    ) = parser.parse_args()  # Now that the options have been defined, instruct the program to parse the command line

    if options.manual == True or len(sys.argv) == 1:
        print """
	Usage:
        calc_waf.py [-M] [-h] {u_file} {u_clim_file} {u_name}
    {v_file} {v_clim_file} {v_name} {zg_file} {zg_clim_file} {zg_name}
    {output_wafx_file} {output_wafy_file} {hPa} {syear} {eyear} {smon} {emon} {sday} {edy} {cyear}

    Options
            -M -> Display this on-line manual page and exit
            -h -> Display a help/usage message and exit

    Description
            Takes as input the u wind, v wind, geopotential height (zg)
        and their climatologies and outputs the x and y
        components of the wave activity flux.

    Assumptions by tnf_xy_onelevel.f90 script (i.e. hard wired elements)
        The WAF Fortran code is hard wired to produce output on a
        global 2.5 by 2.5 deg grid (i.e. 73 lats, 144 lons).
        It is hard wired that the input data is from the 250 hPa level.
        It is assumed that the input data is three dimensional (time, lat, lon).

    Python (cal_WAF.py) script assumptions and capability
        It assumes the global data is passed as input with any degree resolution,
        which will be going to be converted to 2.5x2.5 deg.
        User can choose their desired pressure level from input files.
        User can choose their desired time (range of dates) from the input files,
        it assumes that daily input files are being passed as input but works well for monthly data also.


    Reference
        Takaya, K., Nakamura, H., 2001.
            A formulation of a phase-independent wave-activity flux
        for stationary and migratory quasi-geostrophic eddies
        on zonally-varying basic flow. J. Atmos. Sci., 58, 608-627

    Environment
            Need to load cdat

        Example by (arulalant@gmail.com, arulalan@cas.iitd.ac.in)
        /opt/cdat/bin/cdat calc_waf.py
         era_interim_u_dailyavg_19820101_20171031.nc
         era_interim_u_daily_climatology_1982_2010_cdo.nc
         u
         era_interim_v_dailyavg_19820101_20171031.nc
         era_interim_v_daily_climatology_1982_2010_cdo.nc
         v
         era_interim_z_dailyavg_19820101_20171031.nc
         era_interim_z_daily_climatology_1982_2010_cdo.nc
         z
         WAF_DATA/wafx_200hPa_%s.nc     # %s will be replaced with year
         WAF_DATA/wafy_200hPa_%s.nc     # %s will be replaced with year
         200                            # extract particular pressure level
         1982                           # start year of needed data from input data
         2017                           # end year of needed data from input data
         3                              # start month of needed data from input data
         6                              # end month of needed data from input data
         1                              # start day of start month of needed data from input data
         30                             # end day of end month of needed data from input data
         2010                           # climatology year of climatology input files


    Purpose : Purpose of the extra 8 arguments to fulfil the extracting time and level constraints from the input files
              and create the year wise output wafx, wafy files. For example, in the above example arguments are given to
    	  extract MAMJ months, daily data from 1982-2017 years input and climatology year is 2010 (which was created
    	  by cdo from year 1982-2010).
    Author
            Damien Irving, 12 Oct 2012.
    Updated by :
            Arulalan.T 29 Jan 2018. <arulalant@gmail.com>
            Works on latest UV-CDAT version as on 2018.
            As per fortran script, both raw data and climatology data time length should be same.
            So we are splitting the calc_waf calculation of every year in a loop.
    Bugs
        Please report any problems to: d.irving@student.unimelb.edu.au
	"""
        sys.exit(0)

    else:

        # Repeat the command line arguments #

        print 'Input u wind: ', args[0]
        print 'Input u wind climatology: ', args[1]
        print 'Input v wind: ', args[3]
        print 'Input v wind climatology: ', args[4]
        print 'Input zg: ', args[6]
        print 'Input zg climatology: ', args[7]
        print 'Ouput wafx: ', args[9]
        print 'Ouput wafy: ', args[10]
        print 'lev: ', args[11]
        print 'syear: ', args[12]
        print 'eyear: ', args[13]
        print 'smon: ', args[14]
        print 'emon: ', args[15]
        print 'sday: ', args[16]
        print 'eday: ', args[17]
        print 'cyear: ', args[18]

    fname_u, fname_uclim, vname_u, fname_v, fname_vclim, vname_v, fname_zg, fname_zgclim, vname_zg, fname_wafx_format, fname_wafy_format, hPa, syear, eyear, smon, emon, sday, eday, cyear = args

    for year in range(int(syear), int(eyear) + 1):
        fname_wafx = fname_wafx_format % year
        fname_wafy = fname_wafy_format % year
        print "Progress", year
        main(fname_u, fname_uclim, vname_u, fname_v, fname_vclim, vname_v, fname_zg, fname_zgclim, vname_zg, fname_wafx,
             fname_wafy, hPa, smon, emon, sday, eday, year, cyear)
