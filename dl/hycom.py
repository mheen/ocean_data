from dataclass import from_downloaded as modeldata_from_downloaded
import log as log
from utilities import convert_time_to_datetime,convert_datetime_to_time
from utilities import get_time_indices,get_variable_name
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import os

def opendap_server(output_dir,main_data_url,start_date,end_date,log_file='dl/hycom.log',variables=['u','v'],i_depths=[0]):
    '''Downloads 3-hourly HYCOM data from OpenDAP server and saves data to daily netcdf file.
    
    Input:
    - output_dir [string]: directory where output netcdf files are saved
    - main_data_url [string]: url to relevant HYCOM OpenDAP server
    - start_date [datetime.date]: date to start download
    - end_date [datetime.date]: date to end download
    - log_file (optional) [string]: file in which log messages are saved (default: dl/hycom.log)
    - variables (optional) [list of strings]: list of strings specifying variables to download, these
      should match the variable names in "input/variables.json" (default: ['u','v'])
    - i_depths (optional) [list of integers]: indices of depth layers to download (default: [0], surface)'''
    n_years = end_date.year-start_date.year
    for n_year in range(n_years):
        year = start_date.year+n_year
        input_path = main_data_url+str(year)
        log.info(log_file,'Loading OpenDap data...')
        netcdf = Dataset(input_path)
        times_org = netcdf['time'][:]
        time_units = netcdf['time'].units
        times = convert_time_to_datetime(times_org,time_units)
        for time in times:
            output_path = output_dir+time.strftime('%Y%m%d')+'.nc'
            if os.path.exists(output_path):
                log.info(log_file,f'File {output_path} already exists, skipping download.')
                continue
            # download daily output
            log.info(log_file,f'Downloading {time.strftime("%d-%m-%Y")}')
            i_times = get_time_indices(times,time)
            hycomdata = modeldata_from_downloaded(netcdf,variables,'hycom',i_times=i_times,i_depths=i_depths)           
            hycomdata.write_to_netcdf(output_dir)
        netcdf.close()
