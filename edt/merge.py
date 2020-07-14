from modeldata import ModelData, netcdf_to_dimension, netcdf_to_quantity, from_local_file
from utilities import get_dir, get_ncfiles_in_dir, get_ncfiles_in_time_range
import log
from netCDF4 import Dataset

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
