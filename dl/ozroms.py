from modeldata import from_downloaded as modeldata_from_downloaded
import log as log
from utilities import get_start_and_end_indices
from utilities import convert_time_to_datetime, convert_datetime_to_time
from utilities import get_time_indices, get_variable_name
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import os

def get_lon_index_from_netcdf(input_path,lon_range,lon_description='lon'):
    lon = Dataset(input_path)[lon_description][0,:]
    start_index,end_index = get_start_and_end_indices(lon,lon_range[0],lon_range[-1])
    return (start_index,end_index)

def get_lat_index_from_netcdf(input_path,lat_range,lat_description='lat'):
    lat = Dataset(input_path)[lat_description][:,0]
    start_index,end_index = get_start_and_end_indices(lat,lat_range[0],lat_range[-1])
    return (start_index,end_index)

def opendap_server(output_dir,input_url,lon_range=None,lat_range=None,variables=['u','v'],i_depths=[0],log_file='dl/ozroms.log'):
    if lon_range is not None:
        lon_description = get_variable_name('ozroms','lon')
        i_lon_start,i_lon_end = get_lon_index_from_netcdf(input_url,lon_range,lon_description=lon_description)
        i_lons = slice(i_lon_start,i_lon_end,1)
    else:
        i_lons = None
    if lat_range is not None:
        lat_description = get_variable_name('ozroms','lat')
        i_lat_start,i_lat_end = get_lat_index_from_netcdf(input_url,lat_range,lat_description=lat_description)
        i_lats = slice(i_lat_start,i_lat_end,1)
    else:
        i_lats = None
    netcdf = Dataset(input_url)
    time_description = get_variable_name('ozroms','time')
    times_org = netcdf[time_description][:].filled(fill_value=np.nan)
    times = convert_time_to_datetime(times_org,netcdf[time_description].units)    
    for time in times:
        output_path = output_dir+time.strftime('%Y%m%d')+'.nc'
        if os.path.exists(output_path):
            log.info(log_file,f'File {output_path} already exists, skipping download.')
            continue
        # download daily output
        log.info(log_file,f'Downloading {time.strftime("%d-%m-%Y")}')
        i_times = get_time_indices(times,time)
        ozromsdata = modeldata_from_downloaded(netcdf,variables,'ozroms',i_times=i_times,i_depths=i_depths,i_lons=i_lons,i_lats=i_lats)           
        ozromsdata.write_to_netcdf(output_dir)
    netcdf.close()
