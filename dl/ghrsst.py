from modeldata import from_downloaded as modeldata_from_downloaded
import log as log
from utilities import convert_time_to_datetime, convert_datetime_to_time
from utilities import get_time_indices, get_variable_name
from utilities import get_ncfiles_from_opendap_catalog
from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np
import os

def opendap_server_daily(output_dir, main_data_url, start_date, end_date, log_file='dl/ghrsst.log', variable='sst'):
    '''Downloads daily GHRSST data from NASA OpenDAP server and saves data to netcdf file.

    Input:
    - output_dir [string]: directory where output netcdf files are saved
    - main_data_url [string]: url to relevant NASA GHRSST OpenDAP server
    - start_date [datetime.date]: date to start download
    - end_date [datetime.date]: date to end download
    - log_file (optional) [string]: file in which log messages are saved (default: dl/modis_aqua.log)
    - variable (optional) [string]: string specifying variable to download, this
      should match the variable names in "input/variables.json" (default: 'sst')'''
    model_variable = get_variable_name('ghrsst', variable)
    n_years = end_date.year-start_date.year+1
    for n_year in range(n_years):
        year = start_date.year+n_year
        t0 = datetime(year,1,1)
        if year % 4 == 0:
            n_days = 366
        else:
            n_days = 365
        for day in range(1,n_days+1): # cycle through number of days in year
            date_str = (t0+timedelta(days=day-1)).strftime('%Y%m%d')
            output_path = f'{output_dir}{date_str}.nc'
            if os.path.exists(output_path):
                log.info(log_file, f'File already exists, skipping download: {output_path}')
                continue
            input_dir = main_data_url+f'{year}/{str(day).zfill(3)}/'
            catalog_url = f'{input_dir}contents.html'
            ncfiles = get_ncfiles_from_opendap_catalog(catalog_url)
            for ncfile in ncfiles:
                input_path = input_dir+ncfile
                log.info(log_file,f'Downloading OpenDap data from: {input_path}')
                netcdf = Dataset(input_path)
                ghrsstdata = modeldata_from_downloaded(netcdf,[variable],'ghrsst',depth_value=0)
                ghrsstdata.write_to_netcdf(output_dir)
                netcdf.close()
