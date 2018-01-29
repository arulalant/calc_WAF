# calc_WAF
Wave Activity Flux
------------------

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

	Assumptions (i.e. hard wired elements)
	    The WAF Fortran code is hard wired to produce output on a
	    global 2.5 by 2.5 deg grid (i.e. 73 lats, 144 lons).
	    It is hard wired that the input data is from the 250 hPa level.
	    It is assumed that the input data is three dimensional (time, lat, lon).

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
