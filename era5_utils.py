# era5_utils.py
import xarray as xr
import pandas as pd
import glob
import os




# Location of data
data_loc = '/g/data/rt52/era5/'
data_loc_pressure = data_loc+'pressure-levels/reanalysis/'
data_loc_single = data_loc+'single-levels/reanalysis/'

data_loc_pv = '/g/data/uc16/era5/potential-vorticity/oper/'

# String for level description in filenames
p_str = 'era5_oper_pl'
s_str = 'era5_oper_sfc'
pv_str = 'era5_pv_oper_an'



def open_era5_month(year: int, month: int,
                    day: int = None, hour: int = None,
                    variables=("u","v","t","z")) -> xr.Dataset:
    """
    Open ERA5 monthly NetCDF files for specified variables and optionally
    select a specific time (day/hour).

    Parameters
    ----------
    directory : str
        Path to the folder containing the NetCDF files.
    year : int
        Year of the desired file.
    month : int
        Month of the desired file (1-12).
    day : int, optional
        Day of the month to select.
    hour : int, optional
        Hour of the day to select (0-23).
    variables : tuple of str
        Variables to load (default: ("u","v","t","z")).

    Returns
    -------
    xr.Dataset
        The opened dataset with selected variables, optionally sliced in time.
    """
    files = []
    for var in variables:
 
        # Try to load as a pressure variable 
        file_directory = data_loc_pressure+var+'/'+str(year)+'/'
        pattern = os.path.join(file_directory, f"{var}_{p_str}_{year}{month:02d}*.nc")
        match = glob.glob(pattern)

        if not match:
   
            # Try to load as a single-level variable            
            file_directory = data_loc_single+var+'/'+str(year)+'/'
            pattern = os.path.join(file_directory, f"{var}_{s_str}_{year}{month:02d}*.nc")
            match = glob.glob(pattern)

            if not match:

                # Try to load as a dynamical tropopause variable            
                file_directory = data_loc_pv+var+'/'+str(year)+'/'
                pattern = os.path.join(file_directory, f"{var}_{pv_str}_{year}{month:02d}*.nc")
                match = glob.glob(pattern)
      
            if not match:
                raise FileNotFoundError(f"No {var} file found for {year}-{month:02d}")

        if len(match) > 1:
            raise ValueError(f"Multiple {var} files matched: {match}")
        files.extend(match)

    # Merge all variables into a single dataset
    ds = xr.open_mfdataset(files, combine="by_coords")

    # Time slice if requested
    if day is not None and hour is not None:
        time_sel = pd.Timestamp(year, month, day, hour)
        ds = ds.sel(time=time_sel)

    return ds

