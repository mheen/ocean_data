from modeldata import ModelData, netcdf_to_dimension, netcdf_to_quantity, from_local_file
from modeldata import Dimension
from utilities import get_dir, get_ncfiles_in_dir, get_ncfiles_in_time_range
from utilities import get_variable_name, get_variable_name_reverse
import log
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import os

def merge_variables_into_netcdf(input_dirs,output_dir,timeformat='%Y%m',log_file='edt/merge_netcdfs.log'):
    main_input_dir = input_dirs[0]
    input_dirs = input_dirs[1:]
    nc_files = get_ncfiles_in_dir(main_input_dir)
    for nc_file in nc_files:
        log.info(log_file,f'Loading data {main_input_dir+nc_file}')
        modeldata = from_local_file(main_input_dir+nc_file)
        date0 = modeldata.time.datetime[0].date()
        for input_dir in input_dirs:
            nc_file_add = get_ncfiles_in_time_range(input_dir,date0,date0,timeformat=timeformat)
            if len(nc_file_add) > 1:
                raise ValueError(f'More than one netcdf file found in date range {date0.strftime(timeformat)}')
            netcdf = Dataset(input_dir+nc_file_add[0])        
            variable_names = list(set(netcdf.variables.keys())-set(['time','depth','lat','lon']))
            for variable_name in variable_names:
                log.info(log_file,f'Adding variable {variable_name} to netcdf')
                variable = netcdf_to_quantity(netcdf,variable_name)
                modeldata.fill_variable(variable_name,variable)
            netcdf.close()
        output_path = modeldata.get_output_path(output_dir,filename_format=timeformat)
        log.info(log_file,f'Saving merged netcdf to {output_path}')
        _ = modeldata.write_to_netcdf(output_dir,filename_format=timeformat)

def get_monthly_means_from_local_files(model : str, input_dir : str, output_dir : str,
                                       filename_format = '%Y%m', variables=None,
                                       i_depths=None, i_lats=None, i_lons=None) -> ModelData:
    ncfiles = get_ncfiles_in_dir(input_dir)
    netcdf = Dataset(input_dir+ncfiles[0])
    depth = netcdf_to_dimension(netcdf,get_variable_name(model,'depth'),new_name='depth',i_use=i_depths)
    lon = netcdf_to_dimension(netcdf,get_variable_name(model,'lon'),new_name='lon',i_use=i_lons)    
    lat = netcdf_to_dimension(netcdf,get_variable_name(model,'lat'),new_name='lat',i_use=i_lons)
    netcdf.close()
    for month in np.arange(1,13):
        # initialise modeldata
        time = Dimension('time',np.array([(datetime(1900,month,1)-datetime(1900,1,1)).days]),'days since 1900-01-01')
        modeldata = ModelData(time,depth,lat,lon)
        for variable in variables:
            variable_model = get_variable_name(model,variable)
            n = 0
            for ncfile in ncfiles:
                ncfile_datestr,_ = os.path.splitext(ncfile)
                if datetime.strptime(ncfile_datestr,filename_format).month == month:
                    netcdf = Dataset(input_dir+ncfile)
                    if n == 0:
                        quantity = netcdf_to_quantity(netcdf,variable_model,new_name=variable,i_depths=i_depths)
                    else:
                        next_quantity = netcdf_to_quantity(netcdf,variable_model,new_name=variable,i_depths=i_depths)
                        quantity.values = np.concatenate((quantity.values,next_quantity.values),axis=0)
                    netcdf.close()
                    next_quantity = None
                    n += 1
            quantity.values = np.nanmean(quantity.values,axis=0)
            modeldata.fill_variable(variable,quantity)
            modeldata.write_to_netcdf(output_dir,filename_format)
