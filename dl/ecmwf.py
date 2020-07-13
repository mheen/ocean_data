from modeldata import from_downloaded
from utilities import get_ecmwf_variable_code,get_n_months,add_month_to_timestamp
import log
from ecmwfapi import ECMWFDataServer
from datetime import date,timedelta
from netCDF4 import Dataset
import os

def ecmwfserver(output_dir,variables,start_date,end_date,
                dx=0.75,level='sfc',
                filename_format='%Y%m',log_file='dl/ecmwf.log'):
    parameters = get_ecmwf_parameters_from_variables(variables)
    temp_output_dir = output_dir+'temp/'
    if not os.path.exists(temp_output_dir):
        os.mkdir(temp_output_dir)
    n_months = get_n_months(start_date,end_date)
    for n in range(n_months):
        dl_date = add_month_to_timestamp(start_date,n)
        output_path = output_dir+dl_date.strftime(filename_format)+'.nc'
        temp_output_path = temp_output_dir+dl_date.strftime(filename_format)+'.nc'
        if os.path.exists(output_path):
            log.info(log_file,f'File already exists, skipping download: {output_path}')
            if os.path.exists(temp_output_path):
                log.info(log_file,f'Removing temp file: {temp_output_path}')
                os.remove(temp_output_path)
            continue
        if not os.path.exists(temp_output_path):
            log.info(log_file,f'Downloading temp file: {temp_output_path}')
            dl_era_interim_moda(temp_output_path,dl_date,parameters,dx=dx,level=level)
        netcdf = Dataset(temp_output_path)
        modeldata = from_downloaded(netcdf,variables,'ecmwf',depth_value=-10.)
        log.info(log_file,f'Saving permanent requested file: {output_path}')
        modeldata.write_to_netcdf(output_dir,filename_format=filename_format)
        log.info(log_file,f'Removing temp file: {temp_output_path}')
        os.remove(temp_output_path)

def get_ecmwf_parameters_from_variables(variables):
    parameters = []
    for i,variable in enumerate(variables):
        parameters.append(get_ecmwf_variable_code(variable))
        if i < len(variables)-1:
            parameters.append('/')
    return ''.join(parameters)

def dl_era_interim_moda(temp_output_path,date,parameters,dx=0.75,level='sfc'):    
    server = ECMWFDataServer()
    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": date.strftime('%Y%m%d'),
        "expver": "1",
        "grid": f"{dx}/{dx}",
        "levtype": level,
        "param": parameters,
        "stream": "moda",
        "type": "an",
        "format": "netcdf",
        "target": temp_output_path
    })
