from modeldata import from_downloaded as modeldata_from_downloaded
import log
from utilities import get_ncfiles_in_dir
from utilities import get_urls,get_logins
from utilities import get_n_months,add_month_to_timestamp,convert_time_to_datetime
from netCDF4 import Dataset
import numpy as np
import os
from datetime import datetime
from time import time as time_now
from ftplib import FTP,error_temp

def ftp_server_retry(output_dir,main_ftp_dir,start_date,end_date,temporal_resolution,log_file='dl/cmems.log',variables=['u','v'],i_depths=[0],i_times=None):
    '''Calls the ftp_server function to download CMEMS data from FTP, if the connection to the FTP breaks during
    execution, the function is called again to retry connecting.
    
    Input:
    - output_dir [string]: local directory where final files are stored
    - main_ftp_dir [string]: directory on FTP server where data to download is stored
    - start_date [datetime.date]: date to start download
    - end_date [datetime.date]: date to end download
    - log_file (optional) [string]: file in which log messages are saved (default: dl/cmems.log)
    - variables (optional) [list of strings]: list of strings specifying variables to download, these
      should match the variable names in "input/variables.json" (default: ['u','v'])
    - i_depths (optional) [list of integers]: indices of depth layers to download (default: [0], surface)'''
    try:
        ftp_server(output_dir,main_ftp_dir,start_date,end_date,temporal_resolution,log_file=log_file,variables=variables,i_depths=i_depths)
    except DisconnectedError as e:
        log.error(log_file, 'FTP connection broken: ' + str(e))
        log.info(log_file,'Trying to reconnect to FTP server...')
        ftp_server_retry(output_dir,main_ftp_dir,start_date,end_date,temporal_resolution,log_file=log_file,variables=variables,i_depths=i_depths)

def ftp_server(output_dir,main_ftp_dir,start_date,end_date,temporal_resolution,log_file='dl_cmems.log',variables=['u','v'],i_depths=[0],i_times=None):
    '''Connects to the CMEMS FTP server and downloads full netcdf files from the FTP with either daily or monthly resolution.
    Full netcdf files are stored temporarily on a local drive before requested variables and depth layers are extracted. If this
    is successful, the temporary full netcdf file is then removed. If the connection to the FTP server breaks during execution,
    a DisconnectedError is returned.'''
    ftp_url = get_urls('cmems_main')
    login = get_logins('cmems')
    ftp = connect_to_cmems_ftp(ftp_url,main_ftp_dir,login['user'],login['pwd'],log_file)
    if temporal_resolution == 'monthly':
        _dl_monthly(ftp,output_dir,start_date,end_date,log_file,variables,i_depths,i_times)
    elif temporal_resolution == 'daily':
        _dl_daily(ftp,output_dir,start_date,end_date,log_file,variables,i_depths,i_times)
    else:
        raise ValueError('Unknown temporal resolution requested for CMEMS FTP download. Valid options are "daily" or "monthly".')
    ftp.quit()

def connect_to_cmems_ftp(ftp_url,main_ftp_dir,user,pwd,log_file):
    '''Connects to CMEMS FTP'''
    log.info(log_file,'Accessing CMEMS FTP server: '+ftp_url+'/'+main_ftp_dir)
    ftp = FTP(ftp_url)
    ftp.login(user=user,passwd=pwd)
    ftp.cwd(main_ftp_dir)
    return ftp

def _dl_monthly(ftp,output_dir,start_date,end_date,log_file,variables,i_depths,i_times):
    '''Downloads full monthly netcdf files to temporary directory, saves requested
    variables and depth layers, and removes temporary full netcdf file when successfull.'''
    n_years = end_date.year-start_date.year+1
    # loop over years because netcdf files are stored in year/ directories
    for n_year in range(n_years):
        year = start_date.year+n_year
        ftp.cwd(str(year)+'/')
        log.info(log_file,f'Changed FTP working dir to: {ftp.pwd()}')
        filenames = ftp.nlst()        
        _dl_file(ftp,filenames,output_dir,log_file,variables,i_depths,i_times,filename_format='%Y%m')
        ftp.cwd('../')

def _dl_daily(ftp,output_dir,start_date,end_date,log_file,variables,i_depths,i_times):
    '''Downloads full daily netcdf files to temporary directory, saves requested
    variables and depth layers, and removes temporary full netcdf file when successful.'''
    n_months = get_n_months(start_date,end_date)
    # loop over months because netcdf files are stored in year/month/ directories
    for n_month in range(n_months):
        date = add_month_to_timestamp(start_date,n_month)
        date_dir = str(date.year)+'/'+str(date.month).zfill(2)+'/'
        ftp.cwd(date_dir)
        log.info(log_file,f'Changed FTP working dir to: {ftp.pwd()}')
        filenames = ftp.nlst()        
        _dl_file(ftp,filenames,output_dir,log_file,variables,i_depths,i_times)
        ftp.cwd('../../')

def _dl_file(ftp,filenames,output_dir,log_file,variables,i_depths,i_times,filename_format='%Y%m%d'):
    for filename in filenames:
        # download full file to temporary directory
        temp_output_dir = output_dir+'temp/'
        if not os.path.exists(temp_output_dir):
            os.mkdir(temp_output_dir)
        temp_output_path = temp_output_dir+filename        
        try:
            if not os.path.exists(temp_output_path):
                # skip download of temp file if permanent file already exists
                if filename_format == '%Y%m':
                    output_path = output_dir+filename[-9:]
                elif filename_format == '%Y%m%d':
                    output_path = output_dir+filename[-11:]
                else:
                    raise ValueError('Unknown filename format. Valid options are "%Y%m" and "%Y%m%d".')
                if os.path.exists(output_path):
                    err_message = _check_permanent_netcdf(output_path)
                    if err_message:
                        log.error(log_file,err_message)
                        continue
                    log.info(log_file,f'Permanent file already exists, skipping: {output_path}')
                    continue
                log.info(log_file,f'Downloading temp file: {filename}')
                filehandle = open(temp_output_path,'wb')
                ftp.retrbinary('RETR %s' % filename,filehandle.write)
                filehandle.close()
        except error_temp as e:
            os.remove(temp_output_path)
            raise DisconnectedError("FTP connection broken: " + str(e))
        except Exception as e:
            os.remove(temp_output_path)
            log.error(log_file, str(e), e)
            continue
        if os.path.exists(temp_output_path):
            log.info(log_file, f'File already exists, skipping download: {temp_output_path}')
            # save requested variables and depth levels and remove temporary full file            
            netcdf = Dataset(temp_output_path)
            cmemsdata = modeldata_from_downloaded(netcdf,variables,'cmems',i_times=i_times,i_depths=i_depths)
            output_path = cmemsdata.get_output_path(output_dir,filename_format=filename_format)
            if os.path.exists(output_path):
                err_message = _check_permanent_netcdf(output_path)
                if err_message:
                    log.error(log_file,err_message)
                    continue
                log.info(log_file,f'Permanent file already exists, skipping: {output_path}')
                netcdf.close()
                continue
            log.info(log_file,f'Saving permanent requested file: {output_path}')
            _ = cmemsdata.write_to_netcdf(output_dir,filename_format=filename_format)
            err_message = _check_permanent_netcdf(output_path)
            if err_message:
                log.error(log_file,err_message)
                continue
            log.info(log_file,f'Removing temp file: {temp_output_path}')
            os.remove(temp_output_path)
            os.rmdir(temp_output_dir)
            netcdf.close()

def _check_permanent_netcdf(input_path:str) -> str:    
    # 1. check if netcdf file really exists
    if not os.path.exists(input_path):
        return f'File does not exist: {input_path}'
    # 2. check if netcdf file contains bytes
    if not os.stat(input_path).st_size > 100:
        return f'File contains no bytes: {input_path}'
    # # 3. check if time in netcdf file matches time string in filename
    # filename,_ = os.path.splitext(os.path.basename(input_path))
    # filename_time = datetime.strptime(filename,'%Y%m%d')
    # netcdf = Dataset(input_path)
    # nc_time = convert_time_to_datetime(netcdf['time'],netcdf['time'].units)
    # if not filename_time.date() == nc_time[0].date():
    #     return f'Time in netcdf file does not match name of netcdf file: {input_path}'
    # return None if all checks passed
    return None

class DisconnectedError(RuntimeError):
    def __init__(self, message):
        super().__init__(message)
