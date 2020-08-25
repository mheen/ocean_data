from modeldata import from_downloaded as modeldata_from_downloaded
import log as log
from utilities import convert_time_to_datetime,convert_datetime_to_time
from utilities import get_time_indices,get_variable_name
from utilities import get_ncfiles_from_opendap_catalog
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import os

def opendap_server_monthly(output_dir,main_data_url,start_year,end_year,log_file='dl/modis_aqua.log',variable='sst',horizontal_res=9):
    '''Downloads monthly MODIS-Aqua data from NASA OpenDAP server and saves data to daily netcdf file.
    
    Input:
    - output_dir [string]: directory where output netcdf files are saved
    - main_data_url [string]: url to relevant NASA MODIS-Aqua OpenDAP server
    - start_date [datetime.date]: date to start download
    - end_date [datetime.date]: date to end download
    - log_file (optional) [string]: file in which log messages are saved (default: dl/modis_aqua.log)
    - variable (optional) [string]: string specifying variable to download, this
      should match the variable names in "input/variables.json" (default: 'sst')'''
    model_variable = get_variable_name('modis-aqua',variable)
    n_years = end_year-start_year+1
    for n_year in range(n_years):
        year = start_year+n_year
        for month in range(1,13):
            folder_date = (datetime(year,month,1)).strftime('%m%d')
            input_dir = f'{main_data_url}{year}/{folder_date}/'
            all_ncfiles = get_ncfiles_from_opendap_catalog(input_dir)
            ncfiles = [ncfile for ncfile in all_ncfiles if f'MO.{variable.upper()}.{variable}.{np.int(horizontal_res)}km.nc' in ncfile]
            for ncfile in ncfiles:
                input_path = input_dir+ncfile
                log.info(log_file,f'Downloading OpenDap data from: {input_path}')
                netcdf = Dataset(input_path)
                time_org = [datetime(year,month,1)]
                time,units = convert_datetime_to_time(time_org)
                modisdata = modeldata_from_downloaded(netcdf,[variable],'modis-aqua',depth_value=0,time_values=time,time_units=units)
                modisdata.write_to_netcdf(output_dir)
                netcdf.close()
