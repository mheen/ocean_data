from modeldata import from_downloaded as modeldata_from_downloaded
import log as log
from utilities import get_start_and_end_indices, get_urls, get_dir
from utilities import get_ncfiles_in_dir
from utilities import get_ncfiles_from_opendap_catalog
from utilities import get_n_months, add_month_to_timestamp
from utilities import convert_time_to_datetime, convert_datetime_to_time
from utilities import get_time_indices, get_variable_name
from netCDF4 import Dataset
from datetime import datetime, timedelta
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

def get_ozroms_daily_ncfiles_from_irds(irds_mount,start_date,end_date):
    all_ncfiles = get_ncfiles_in_dir(irds_mount)
    ncfiles = get_ozroms_daily_ncfiles(all_ncfiles,start_date,end_date)
    return ncfiles

def get_ozroms_monthly_ncfiles_from_opendap(catalog_url,start_date,end_date):
    all_ncfiles = get_ncfiles_from_opendap_catalog(catalog_url)
    ncfiles = get_ozroms_monthly_ncfiles(all_ncfiles,start_date,end_date)
    return ncfiles

def get_ozroms_monthly_ncfiles(all_ncfiles,start_date,end_date):
    ncfiles = []
    n_months = get_n_months(start_date,end_date)
    for n in range(n_months):
        time = add_month_to_timestamp(start_date,n)
        for ncfile in all_ncfiles:
            if ncfile[:-3].endswith(str(time.year)):
                month = time.strftime('%b') # 3 letter month
                if ncfile[:-8].endswith(month):
                    ncfiles.append(ncfile)
    return ncfiles

def get_ozroms_daily_ncfiles(all_ncfiles,start_date,end_date):
    ncfiles = []
    n_days = (end_date-start_date).days
    for n in range(n_days):
        time = start_date+timedelta(days=n)
        for ncfile in all_ncfiles:
            if ncfile[:-3].endswith(str(time.year)):
                month = time.strftime('%b')
                day = time.strftime('%d')
                if ncfile[:13].endswith(f'{month}_{day}'):
                    ncfiles.append(ncfile)
    return ncfiles

def _get_i_lons_and_i_lats(input_url,lon_range,lat_range):
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
    return i_lons,i_lats

def bathymetry(input_path, output_path, lon_range=None, lat_range=None):
    i_lons,i_lats = _get_i_lons_and_i_lats(input_path, lon_range, lat_range)
    netcdf = Dataset(input_path)
    modeldata = modeldata_from_downloaded(netcdf, ['h'], 'ozroms', i_lats=i_lats, i_lons=i_lons)
    modeldata.write_to_netcdf(None, output_path=output_path)

def hourly_date_range_from_irds_file(input_path, output_dir, date_range, lon_range=None, lat_range=None,
                                     variables = ['u', 'v'], i_depths=[0], log_file='dl/ozroms_hourly.log'):
    log.info(log_file, f'Accessing: {input_path}')
    i_lons,i_lats = _get_i_lons_and_i_lats(input_path, lon_range, lat_range)
    netcdf = Dataset(input_path)
    time_description = get_variable_name('ozroms', 'time')
    times_org = netcdf[time_description][:].filled(fill_value=np.nan)
    times = convert_time_to_datetime(times_org, netcdf[time_description].units)
    n_days = (date_range[1]-date_range[0]).days
    for i in range(n_days):
        time = date_range[0]+timedelta(days=i)
        output_path = f'{output_dir}{time.strftime("%Y%m%d")}.nc'
        if os.path.exists(output_path):
            log.info(log_file, f'File already exists, skipping: {output_path}')
            continue
        log.info(log_file, f'Downloading {time.strftime("%d-%m-%Y")}')
        i_times = get_time_indices(times, time)
        ozromsdata = modeldata_from_downloaded(netcdf, variables, 'ozroms',
                                               i_times=i_times, i_depths=i_depths,
                                               i_lons=i_lons, i_lats=i_lats)
        ozromsdata.write_to_netcdf(output_dir)
    netcdf.close()

def daily_date_range_from_irds_server(output_dir,date_range,lon_range=None,lat_range=None,
                                variables=['u','v'],i_depths=[0],
                                irds_mount=get_dir('ozroms_daily_input'),log_file='dl/ozroms_daily.log'):
    log.info(log_file,f'Accessing IRDS: {irds_mount}')
    ncfiles = get_ozroms_daily_ncfiles_from_irds(irds_mount,date_range[0],date_range[1])
    log.info(log_file,f'Found daily ncfiles: {ncfiles[0:5]}')
    for ncfile in ncfiles:
        input_path = f'{irds_mount}{ncfile}'
        file_from_opendap_server(output_dir,input_path,lon_range=lon_range,
                                 lat_range=lat_range,variables=variables,i_depths=i_depths,
                                 log_file=log_file)

def date_range_from_opendap_server(output_dir,date_range,lon_range=None,lat_range=None,
                                   variables=['u','v'],i_depths=[0],
                                   catalog_url=get_urls('ozroms_catalog'),
                                   main_url=get_urls('ozroms_main'),
                                   log_file='dl/ozroms_daily.log'):
    ncfiles = get_ozroms_monthly_ncfiles_from_opendap(catalog_url,date_range[0],date_range[-1])
    for ncfile in ncfiles:
        input_url = f'{main_url}{ncfile}'
        file_from_opendap_server(output_dir,input_url,lon_range=lon_range,
                                 lat_range=lat_range,variables=variables,i_depths=i_depths,
                                 log_file=log_file)

def file_from_opendap_server(output_dir,input_url,lon_range=None,lat_range=None,
                             variables=['u','v'],i_depths=[0],log_file='dl/ozroms.log'):
    i_lons,i_lats = _get_i_lons_and_i_lats(input_url,lon_range,lat_range)
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
