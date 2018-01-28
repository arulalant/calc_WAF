# calc_WAF
Wave Activity Flux
------------------

    usage = "usage: %prog [options] {u_file} {u_name} {v_file} {v_name} {output_file}"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-M", "--manual",action="store_true",dest="manual",default=False,help="output a detailed description of the program")

    (options, args) = parser.parse_args()            # Now that the options have been defined, instruct the program to parse the command line

    if options.manual == True or len(sys.argv) == 1:
	print """
	Usage:
            calc_waf.py [-M] [-h] {u_file} {u_clim_file} {u_name} 
	    {v_file} {v_clim_file} {v_name} {zg_file} {zg_clim_file} {zg_name}
	    {output_wafx_file} {output_wafy_file}

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

        Example (abyss.earthsci.unimelb.edu.au)
	    /opt/cdat/bin/cdat calc_waf.py 
	    /work/dbirving/datasets/Merra/data/ua_Merra_250hPa_monthly_native.nc
	    /work/dbirving/datasets/Merra/data/processed/ua_Merra_250hPa_monthly-clim-1981-2010_native.nc
	    ua
	    /work/dbirving/datasets/Merra/data/va_Merra_250hPa_monthly_native.nc
	    /work/dbirving/datasets/Merra/data/processed/va_Merra_250hPa_monthly-clim-1981-2010_native.nc
	    va
	    /work/dbirving/datasets/Merra/data/zg_Merra_250hPa_monthly_native.nc
	    /work/dbirving/datasets/Merra/data/processed/zg_Merra_250hPa_monthly-clim-1981-2010_native.nc
	    zg
            /work/dbirving/datasets/Merra/data/processed/wafx_Merra_250hPa_monthly_native.nc
	    /work/dbirving/datasets/Merra/data/processed/wafy_Merra_250hPa_monthly_native.nc
	    
	Author
            Damien Irving, 12 Oct 2012.

	Bugs 
	    Please report any problems to: d.irving@student.unimelb.edu.au
	"""
