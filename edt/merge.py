from modeldata import ModelData, netcdf_to_dimension, netcdf_to_quantity, from_local_file
from modeldata import Dimension
from utilities import get_dir, get_ncfiles_in_dir, get_ncfiles_in_time_range
from utilities import get_variable_name, get_variable_name_reverse
import log
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import os

def merge_netcdfs_in_folder_to_single_file(input_dir,output_path,variables=['u','v'],velocities=None,log_file=None):
    ncfiles = get_ncfiles_in_dir(input_dir)
    log.info(log_file,f'Loading data {input_dir}{ncfiles[0]}')
    modeldata = from_local_file(input_dir+ncfiles[0],variables=variables)
    for ncfile in ncfiles[1:]:
        netcdf = Dataset(input_dir+ncfile)        
        for variable_name in variables:
            log.info(log_file,f'Appending {variable_name} data from {input_dir+ncfile}')
            modeldata.append_to_variable(variable_name,netcdf[variable_name][:].filled(fill_value=np.nan))
        log.info(log_file,f'Appending time data from {input_dir+ncfile}')
        modeldata.append_to_time(netcdf['time'][:].filled(fill_value=np.nan))
        netcdf.close()
    if velocities == 'currents':
        log.info(log_file,f'Adding absolute current velocity')
        modeldata.add_absolute_current_velocity()
    elif velocities == 'wind':
        log.info(log_file,f'Adding absolute wind velocity')
        modeldata.add_absolute_wind_velocity()    
    modeldata.write_to_netcdf(None,output_path=output_path)

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
            if len(nc_file_add) == 0:
                continue
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

def get_monthly_means_from_local_daily_files(input_dir : str, output_dir : str,
                                             input_filenameformat = '%Y%m%d',
                                             output_filenameformat = '%Y%m',
                                             log_file='merge_daily_to_monthly.log'):
    ncfiles = get_ncfiles_in_dir(input_dir)
    for ncfile in ncfiles:
        ncfile_datestr,_ = os.path.splitext(ncfile)
        ncfile_datetime = datetime.strptime(ncfile_datestr,input_filenameformat)
        output_path = output_dir+ncfile_datetime.strftime(output_filenameformat)+'.nc'
        if os.path.exists(output_path):
            log.info(log_file,f'File already exists, skipping: {output_path}')
            continue
        log.info(log_file,f'Getting monthly means for {ncfile_datetime.month} {ncfile_datetime.year}')
        start_date = datetime(ncfile_datetime.year,ncfile_datetime.month,1)
        if ncfile_datetime.month < 12:
            end_date = datetime(ncfile_datetime.year,ncfile_datetime.month+1,1)
            ncfiles_merge = get_ncfiles_in_time_range(input_dir,start_date,end_date,including_end=0)
        else:
            end_date = datetime(ncfile_datetime.year,ncfile_datetime.month,31)
            ncfiles_merge = get_ncfiles_in_time_range(input_dir,start_date,end_date,including_end=1)
        log.info(log_file,f'Loading data from {input_dir+ncfiles_merge[0]}')
        modeldata = from_local_file(input_dir+ncfiles_merge[0])
        for i in range(1,len(ncfiles_merge)):
            netcdf = Dataset(input_dir+ncfiles_merge[i])
            variables = list(set(netcdf.variables.keys())-set(['time','depth','lat','lon']))
            for variable in variables:
                log.info(log_file,f'Appending data from {input_dir+ncfiles_merge[i]}')
                variable_values = netcdf[variable][:].filled(fill_value=np.nan)
                modeldata.append_to_variable(variable,variable_values)
            netcdf.close()
        log.info(log_file,f'Getting time mean values')
        modeldata.take_time_mean()
        modeldata.get_output_path(output_dir,filename_format=output_filenameformat)
        log.info(log_file,f'Writing to {output_path}')
        modeldata.write_to_netcdf(output_dir,filename_format=output_filenameformat)

def get_total_mean_from_local_files(input_dir,output_path,i_depths=None,i_lats=None,i_lons=None):
    ncfiles = get_ncfiles_in_dir(input_dir)
    netcdf = Dataset(input_dir+ncfiles[0])
    time = Dimension('time',np.array([0.0]),'days since 1900-01-01')
    depth = netcdf_to_dimension(netcdf,'depth',new_name='depth',i_use=i_depths)
    lon = netcdf_to_dimension(netcdf,'lon',new_name='lon',i_use=i_lons)    
    lat = netcdf_to_dimension(netcdf,'lat',new_name='lat',i_use=i_lons)
    variables = list(set(netcdf.variables.keys())-set(['time','depth','lat','lon']))
    netcdf.close()
    modeldata = ModelData(time,depth,lat,lon)
    for variable in variables:
        n = 0
        for ncfile in ncfiles:
            netcdf = Dataset(input_dir+ncfile)
            if n == 0:
                quantity = netcdf_to_quantity(netcdf,variable,new_name=variable,i_depths=i_depths)
            else:
                next_quantity = netcdf_to_quantity(netcdf,variable,new_name=variable,i_depths=i_depths)
                quantity.values = np.concatenate((quantity.values,next_quantity.values),axis=0)
            netcdf.close()
            next_quantity = None
            n += 1
        quantity.values = np.nanmean(quantity.values,axis=0)
        modeldata.fill_variable(variable,quantity)
        modeldata.write_to_netcdf(None,output_path=output_path)

def get_total_monthly_mean_from_local_files(model : str, input_dir : str, output_dir : str,
                                       filename_format = '%Y%m', variables=None,
                                       i_depths=None, i_lats=None, i_lons=None) -> ModelData:
    ncfiles = get_ncfiles_in_dir(input_dir)
    netcdf = Dataset(input_dir+ncfiles[0])
    depth = netcdf_to_dimension(netcdf,get_variable_name(model,'depth'),new_name='depth',i_use=i_depths)
    lon = netcdf_to_dimension(netcdf,get_variable_name(model,'lon'),new_name='lon',i_use=i_lons)    
    lat = netcdf_to_dimension(netcdf,get_variable_name(model,'lat'),new_name='lat',i_use=i_lons)
    if variables is None:
        variables = list(set(netcdf.variables.keys())-set(['time','depth','lat','lon']))
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
