from modeldata import from_downloaded as modeldata_from_downloaded
import log as log
from utilities import convert_time_to_datetime, convert_datetime_to_time
from utilities import get_time_indices, get_variable_name
from utilities import get_ncfiles_from_opendap_catalog
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import os


def opendap_server_monthly(output_dir, main_data_url, start_year, end_year, log_file='dl/modis_aqua.log', variable='sst', horizontal_res=9):
    '''Downloads monthly MODIS-Aqua data from NASA OpenDAP server and saves data to daily netcdf file.

    Input:
    - output_dir [string]: directory where output netcdf files are saved
    - main_data_url [string]: url to relevant NASA MODIS-Aqua OpenDAP server
    - start_date [datetime.date]: date to start download
    - end_date [datetime.date]: date to end download
    - log_file (optional) [string]: file in which log messages are saved (default: dl/modis_aqua.log)
    - variable (optional) [string]: string specifying variable to download, this
      should match the variable names in "input/variables.json" (default: 'sst')'''
    model_variable = get_variable_name('modis-aqua', variable)
    n_years = end_year-start_year+1
    for n_year in range(n_years):
        year = start_year+n_year
        for month in range(1, 13):
            output_path = f'{output_dir}{datetime(year,month,1).strftime("%Y%m%d")}.nc'
            if os.path.exists(output_path):
                log.info(
                    log_file, f'File already exists, skipping download: {output_path}')
                continue
            if variable == 'sst':
                folder_date = (datetime(year, month, 1)).strftime('%m%d')
                input_dir = f'{main_data_url}{year}/{folder_date}/'
            else:
                day_of_year = datetime(year, month, 1).timetuple().tm_yday
                input_dir = f'{main_data_url}{year}/{str(day_of_year).zfill(3)}/'
            all_ncfiles = get_ncfiles_from_opendap_catalog(input_dir)
            if variable == 'sst':
                variable_string = f'MO.SST.sst.{np.int(horizontal_res)}km.nc'
            elif variable == 'chl_a':
                variable_string = f'MO_CHL_chlor_a_{np.int(horizontal_res)}km.nc'
            elif variable == 'turbidity':
                variable_string = f'MO_KD490_Kd_490_{np.int(horizontal_res)}km.nc'
            else:
                raise ValueError(
                    'Unknown variable requested for MODIS-Aqua download. Valid options are: sst, chl_a, turbidity.')
            ncfiles = [
                ncfile for ncfile in all_ncfiles if variable_string in ncfile]
            for ncfile in ncfiles:
                input_path = input_dir+ncfile
                log.info(
                    log_file, f'Downloading OpenDap data from: {input_path}')
                netcdf = Dataset(input_path)
                time_org = [datetime(year, month, 1)]
                time, units = convert_datetime_to_time(time_org)
                modisdata = modeldata_from_downloaded(
                    netcdf, [variable], 'modis-aqua', depth_value=0, time_values=time, time_units=units)
                modisdata.write_to_netcdf(output_dir)
                netcdf.close()
