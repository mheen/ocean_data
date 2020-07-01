from dataclass import DataConverter
import log
from utilities import get_ncfiles_in_dir,get_variable_name,get_variable_name_reverse
from utilities import convert_time_to_datetime,get_n_months,get_l_time_range,add_month_to_timestamp
from netCDF4 import Dataset
from datetime import datetime,timedelta
import numpy as np
import os

def variables_at_depth(input_dir,output_dir,model_name,i_times=[0],i_depths=[0],variables=['u','v','mld','salinity','temp','sea_ice_cover'],filename_format='%Y%m%d',log_file='edt/extract_depth.log'):    
    filenames = get_ncfiles_in_dir(input_dir)
    for filename in filenames:
        data = Dataset(input_dir+filename)
        extracted_data = DataConverter(data,i_times,i_depths,variables,model_name)
        output_path = extracted_data.get_output_path(output_dir,filename_format=filename_format)
        if os.path.exists(output_path):
            log.info(log_file,f'File already exists, skipping: {output_path}')
            data.close()
            continue
        log.info(log_file,f'Extracting depth layers {str(i_depths)} and saving to file: {output_path}')
        _ = extracted_data.write_to_netcdf(output_dir,filename_format=filename_format)
        data.close()

def monthly_files(input_path,output_dir,model_name,log_file='edt/extract_monthly.log'):
    data = Dataset(input_path)
    i_depths = np.arange(len(data['depth'][:])) # keep original depth layers
    # keep original values: find general name from model specific name
    model_variables = list(set(data.variables.keys())-set(data.dimensions.keys()))
    variables = []
    for model_variable in model_variables:
        variables.append(get_variable_name_reverse(model_name,model_variable))
    time_var = get_variable_name(model_name,'time')
    time = convert_time_to_datetime(data[time_var][:],data[time_var].units)
    n_months = get_n_months(time[0],time[-1])
    date0 = datetime(time[0].year,time[0].month,1)
    for i in range(n_months):
        start_date = add_month_to_timestamp(date0,i)
        end_date = add_month_to_timestamp(date0,i+1)-timedelta(days=1)
        l_times = get_l_time_range(time,start_date,end_date)
        i_times = np.where(l_times)[0]
        extracted_data = DataConverter(data,i_times,i_depths,variables,model_name)
        output_path = extracted_data.get_output_path(output_dir,filename_format='%Y%m')
        if os.path.exists(output_path):
            log.info(log_file,f'File already exists, skipping: {output_path}')
            continue
        log.info(log_file,f'Writing monthly data to file: {output_path}')
        _ = extracted_data.write_to_netcdf(output_dir,filename_format='%Y%m')
    data.close()
